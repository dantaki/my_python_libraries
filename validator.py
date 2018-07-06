#!/isr/env python
import pysam,sys,argparse,os,re
from pybedtools import BedTool as Bed
import numpy as np
def argP():
	parser = argparse.ArgumentParser(description="validate SV with coverage and supporting reads",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-i',help='bam',required=True,type=str)
	parser.add_argument('-r', help='restrict search to region',type=str,default=None)
	parser.add_argument('-R', help='restrict to regions listed in the file',type=str,default=None)
	parser.add_argument('-C', help='median chromosome coverage file', type=str, default=None, required=True) 
	parser.add_argument('-v', help='print out SV supporting read names',default=False,action="store_true",required=False)
	parser.add_argument('-o', help='outputfile',required=False,type=str,default="validator.out")
	return parser.parse_args()
def fileErr(fh):
	if not os.path.isfile(fh):
		sys.stderr.write('ERROR: {} CANNOT BE FOUND\n'.format(fh))
		sys.exit(1)
def loadCov(i):
	cov={}
	with open(i) as f:
		for l in f:
			iid, chrom,covv = l.rstrip('\n').split('\t')
			chrom= chrom.replace('chr','')
			cov[(iid,chrom)]=float(covv)
	return cov
def getIID(i):
	iid = i.split('/').pop()
	iid = iid.replace('_rmdup_sorted_rgged_bqsr.bam','')
	iid = iid.replace('-sorted-rmdups-realigned-bqsr-NoBIBD.bam','')
	iid = iid.replace('.bam','')
	if '_' in iid:
		fid, iid = iid.split('_',1)
	return iid
def depth(bamfh,region,cut):
	depth=[]
	c,s,e = splitRegion(region)
	for c,s,e in Bed('chr{} {} {}'.format(c,s,e),from_string=True).subtract("/home/dantakli/sperm/plots/hg19_excluded.bed.gz"):
		reg = '{}:{}-{}'.format(c,s,e)
		reg = reg.replace('chr','')
		for x in realdepth(bamfh,reg,cut): depth.append(x)
	print (region,'DEPTH:',np.median(depth),np.mean(depth))
	return np.median(depth)	
def realdepth(bamfh,region,cut):
	depth=[]
	#depth_result = pysam.depth("-a","-r",region,bamfh)
	depth_result = pysam.depth("-a","-Q" "40","-r",region,bamfh)
	str_flag=0
	if isinstance(depth_result,str):
		depth_result = depth_result.split('\n')
		str_flag=1
	for x in depth_result:
		r = x.rstrip('\n').split('\t')
		if str_flag == 1:
			if len(r)!=3: continue
			depth.append(float(r[2]))	
		else: depth.append(float(r[2]))
	return depth
def discordant(bam,region,sz):
	READS=[]
	if sz < 0: return READS
	else:
		for al in bam.fetch(region=region):
			if 0.925 <= abs(al.template_length)/float(sz) <= 1.025: 
				_left = al.reference_start
				_right = al.reference_end
				if _left == None or _right == None: continue
				ref_pos = list(xrange(_left,_right))
				med_ind = int(len(ref_pos)/2.)
				READS.append( [al.query_name,ref_pos[med_ind] ])
		return READS
class Cigar():
	def __init__(self,cig=None):
		aLen = 0
		"""----------------------------------------------------------------"""
		"""Courtesy of Pysam:http://pysam.readthedocs.io/en/latest/api.html"""
		CODE2CIGAR= ['M','I','D','N','S','H','P','=','X','B']
		#if PY_MAJOR_VERSION >= 3:
			#CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
		#else:
		CIGAR2CODE = dict([ord(y), x] for x, y in enumerate(CODE2CIGAR))
		CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
		parts = CIGAR_REGEX.findall(cig)
		self.cig=[(CIGAR2CODE[ord(y)], int(x)) for x,y in parts]
		"""----------------------------------------------------------------"""
		left = 0
		right = 0
		leftFlg, leftLen = self.cig[0]
		rightFlg, rightLen = self.cig[-1]
		if leftFlg == 4 or leftFlg == 5: left = leftLen
		if rightFlg == 4 or rightFlg == 5: right = rightLen
		self.leftClip = left
		self.rightClip = right
	def aLen(self):
		aLen=0
		for (flg,leng) in  self.cig:
			if flg == 0 or flg == 2 or flg==3 or flg==7 or flg==8: aLen=aLen+leng
		return aLen
def splitRegion(region):
	chrom, pos = region.split(':')
	s,e = pos.split('-')
	return chrom,s,e
def getBreak(s1,e1,s2,e2,svType):
	a = sorted([s1,e1,s2,e2])
	if svType=='DEL': return a[1],a[2]
	elif svType=='DUP': return a[0],a[3]
def split(bam,region,sv,svType):
	READS=[]
	chrom,s,e = splitRegion(sv)
	for al in bam.fetch(region=region):

		pName,pChrom,pStart,pEnd,pStrand = al.query_name, al.reference_name, al.reference_start, al.reference_end, '+'
		if not al.has_tag('SA'): continue
		if al.is_reverse: pStrand = '-'
		salns = str(al.get_tag('SA')).split(';')
		if not salns[-1].startswith('SA'): del salns[-1]
		_left = al.reference_start
		_right = al.reference_end

		if _left == None or _right == None: continue

		ref_pos = list(xrange(_left,_right))
		med_ind = int(len(ref_pos)/2.)

		for saln in salns:
			salnList = saln.split(',')
			sRef,sLeft,sStrand,sCigar,sQ = salnList[0:5]
			sCig = Cigar(sCigar)
			sLeft=int(sLeft)-1
			sRight = int(sLeft)+sCig.aLen()
			if pChrom == sRef and pStrand == sStrand:
				brkS,brkE = getBreak(pStart,pEnd,sLeft,sRight,svType)
				if Bed('{} {} {}'.format(chrom,s,e),from_string=True).intersect(Bed('{} {} {}'.format(chrom,brkS,brkE),from_string=True),f=0.95,F=0.95).count() > 0: READS.append( [pName,ref_pos[med_ind] ] )
	return READS 
def clipped(bam,leftBreak,rightBreak,sv,svType):
	READS=[]
	chrom,s,e = splitRegion(sv)
	for al in bam.fetch(region=leftBreak):
		pName,pChrom,pStart,pEnd = al.query_name, al.reference_name, al.reference_start, al.reference_end
		if al.cigarstring == None: continue
		cig = Cigar(al.cigarstring)
		if svType == 'DEL':
			if cig.rightClip > 0 and abs(pEnd-int(s)) < 26: READS.append(pName)
		if svType == 'DUP':
			if cig.leftClip > 0 and abs(pStart - int(s)) < 26: READS.append(pName)
	for al in bam.fetch(region=rightBreak):
		pName,pChrom,pStart,pEnd = al.query_name, al.reference_name, al.reference_start, al.reference_end
		if al.cigarstring==None: continue
		cig = Cigar(al.cigarstring)
		if svType == 'DEL':
			if cig.leftClip > 0 and abs(pStart-int(e)) < 26: READS.append(pName)
		if svType == 'DUP':
			if cig.rightClip > 0 and abs(pEnd - int(e)) < 26: READS.append(pName)
	return READS
def fetchRegion(ifh,region,svType,cov):
	"""get depth, get split reads, get clipped reads, get discordant reads"""
	supporting=[]
	reads=[]
	chrom,s,e = splitRegion(region)
	sz = int(e)-int(s)
	flank = 150
	leftBreak = '{}:{}-{}'.format(chrom,int(s)-flank,int(s)+flank)
	rightBreak = '{}:{}-{}'.format(chrom,int(e)-flank,int(e)+flank)
	bam, iid  = pysam.AlignmentFile(ifh,'rb'), getIID(ifh)
	cov_chrom = chrom
	#if 'REACH' not in iid and str(chrom)=='22': cov_chrom='genome'
	sys.stderr.write('IID: {} {}\n'.format(iid,cov[(iid,cov_chrom)]))
 	
	cut = 150
	if 'blood' in iid or 'sperm' in iid: cut = 250

	CN = 2 * depth(ifh,region,cut)/cov[(iid,cov_chrom)]
	for al in bam.fetch(region=leftBreak): reads.append(al.query_name)
	for al in bam.fetch(region=rightBreak): reads.append(al.query_name)
	_left,_right=[],[]
	for x in discordant(bam,leftBreak,sz): 
		supporting.append(x[0])
		_left.append(x[1])
	for x in discordant(bam,rightBreak,sz): 
		supporting.append(x[0])
		_right.append(x[1])
	for x in split(bam,leftBreak,region,svType): 
		supporting.append(x[0])
		_left.append(x[1])
	for x in split(bam,rightBreak,region,svType): 
		supporting.append(x[0])
		_right.append(x[1])

	left_med_brk= s
	right_med_brk = e
	if len(_left) > 0:
		left_med_brk = int(np.median(_left))
	
	if len(_right) >0:
		right_med_brk = int(np.median(_right))

	all_reads= []
	for al in bam.fetch(region='{}:{}-{}'.format(chrom,left_med_brk,left_med_brk)): 
		if al.query_name not in supporting: all_reads.append(al.query_name)
	for al in bam.fetch(region='{}:{}-{}'.format(chrom,right_med_brk,right_med_brk)): 
		if al.query_name not in supporting: all_reads.append(al.query_name)

	#for x in clipped(bam,leftBreak,rightBreak,region,svType): supporting.append(x)
	supporting = list(set(sorted(supporting)))
	reads = list(set(sorted(all_reads)))
	return len(reads),len(supporting),len(supporting)/float(len(reads)+len(supporting)), CN, supporting
	#print region,svType,len(supporting),len(reads),len(supporting)/float(len(reads))
if __name__ == '__main__':
	args = argP()
	ifh = args.i
	region = args.r
	regionFile = args.R
	fileErr(args.C)
	covs = loadCov(args.C)
	ofh = args.o
	fileErr(ifh)
	orn=""
	if args.v == True:orn = open(ofh.replace(".","_readNames."),'w')
	out = open(ofh,'w')
	out.write('IID\tCHROM\tSTART\tEND\tTYPE\tCN\tCONC\tDISC\tPROP\n')
	if regionFile == None: 
		conc, disc, prop, cn = fetchRegion(ifh,region,covs)	
	else:  	
		fileErr(regionFile)
		with open(regionFile,'r') as f:
			for l in f:
				a = l.rstrip('\n').split('\t')
				region = str(a[0]+':'+str(int(a[1])-1)+'-'+a[2])	
				conc, disc, prop, cn, supporting = fetchRegion(ifh,region,a[3],covs)
				out.write('\t'.join(map(str,(getIID(ifh),a[0],a[1],a[2],a[3],cn,conc,disc,prop)))+'\n')
				if args.v == True:
					if len(supporting)>0: orn.write('\n'.join(supporting)+'\n')
					orn.close()
	out.close()


