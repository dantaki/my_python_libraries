import argparse,pysam, os,sys
from khatam import Cigar
# foreach breakpoint
# get split-reads
   # match split-reads to SV
   # check if soft-clipped
   # get MAPQ and sequence
class Args():
	def __init__(self):
		parser = argparse.ArgumentParser(description="Fetch breakpoint sequence",formatter_class=argparse.RawTextHelpFormatter)
		parser.add_argument('-i',help='bam file',required=True,type=str)
		parser.add_argument('-b', help='bed file',type=str,default=None,required=True)
		parser.add_argument('-n',help='number of flanking base pairs to search for split reads [default: 150]',required=False,type=int,default=150)
		parser.add_argument('-x',help='minimum distance to breakpoints [default: 20]',required=False,type=int,default=20)
		parser.add_argument('-o', help='output file',required=False,type=str,default="breaks.out")
		args = parser.parse_args()
		self.bam,self.bed,self.out = args.i,args.b,args.o
		self.flank,self.min_dist = args.n,args.x
class Read():
	def  __init__(self,aln):
		self.qname = aln.query_name
		self.chrom=aln.reference_name
		self.left=aln.reference_start
		self.right=aln.reference_end
		self.cigar = Cigar.Cigar(aln.cigarstring)
		self.cigarstring = aln.cigarstring
		self.strand = '+'
		if aln.is_reverse==True: self.strand='-'
		self.seq = aln.query_alignment_sequence
		self.left_hard,self.right_hard=False,False
		if self.cigar.cig[0][0] == 5: self.left_hard=True
		if self.cigar.cig[-1][0] == 5: self.right_hard=True
		self.mapq = aln.mapping_quality
def err_fh(fh):
	if not os.path.isfile(fh): 
		sys.stderr.write('FATAL ERROR: {} NOT FOUND\n'.format(fh))
		sys.exit(1)

def fetch_reads(bam,regions):
	reads={}
	for region in regions:
		for aln in bam.fetch(region=region):
			left_aln,right_aln = None,None #left and right alignments for the breakpoint 
			if not aln.has_tag('SA'):continue
			read = Read(aln)
			if reads.get(read.qname)==None: reads[read.qname]=[read]
			else: reads[read.qname].append(read)
	return reads
def deletion_breakpoint(reads):
	for read_name in reads:
		left_read, right_read = None, None
		leftmost=None
		if len(reads[read_name])!=2: continue 
		for Read in reads[read_name]:
			if leftmost==None: 
				leftmost=Read.left
				left_read=Read
			else:
				if Read.left < leftmost: 
					leftmost=Read.left
					right_read = left_read
					left_read=Read
				else: 
					right_read=Read
		print left_read.right,right_read.left, left_read.seq, right_read.seq 
		#print read.qname,read.mapq,read.left,read.right,read.left_hard,read.right_hard,read.cigarstring
		#parse_supplementary_alignments(aln,'DEL',start,end,min_dist)
def parse_supplementary_alignments(aln,svtype,start,end,min_dist):
	"""
	check each supplementary alignment for SV evidence
	"""
	prm_strand = '+' # forward
	if aln.is_reverse: prm_strand='-' #reverse
	sec_alignments = str(aln.get_tag('SA')).split(';')
	if not sec_alignments[-1].startswith('SA'): del sec_alignments[-1] # removes empty entry 
	n_aln = len(sec_alignments)+1
	for sec_aln in sec_alignments:
		entry = sec_aln.split(',')
		sec_contig, sec_left_pos, sec_strand, sec_cigar, sec_mapq = entry[0:5]
		# skip if alignments are on not the same strand
		if prm_strand != sec_strand: continue
		# skip if alignments are not on the same chromosome 
		if sec_contig != aln.reference_name:continue
		if svtype=='DEL': 
			left_left,left_right,right_left,right_right = deletion_breakpoints(aln,entry)
			left_dist = abs(left_right-start)
			right_dist = abs(right_left-end)
			if left_dist > min_dist or right_dist > min_dist: continue
			aln_cig = Cigar.Cigar(aln.cigarstring)
			ref_pos = aln.get_reference_positions()
			left_seq, right_seq = None,None
			if aln_cig.right_clip > 0:
				if 'H' not in aln.cigarstring: # Soft clipped
					soft_end = len(aln.query_sequence)-aln_cig.right_clip
					left_seq = aln.query_sequence[0 : soft_end]
					right_seq = aln.query_sequence[soft_end:-1]
					print aln.query_name,n_aln,aln.mapping_quality,aln.cigarstring,left_seq,sec_mapq,sec_cigar,right_seq
			#if aln_cig.right_clip > 0:
			#	left_seq = aln.query_sequence[aln.query_alignment_start:aln.query_alignment_end]
			#	right_seq = aln.query_sequence[aln.query_alignment_end-1:-1]
				#print left_seq,'     ',right_seq
			#print aln.query_name, aln.query_sequence, aln.cigarstring,prm_strand,aln.query_alignment_sequence

def deletion_breakpoints(aln,entry):
	"""
	check if the two alignments support a deletion

	input: primary alignment object from AlignmentFile().fetch(), list of a SA tag entry
	"""
	left_cigar,left_left,right_cigar,right_left = aln.cigarstring,aln.reference_start,entry[3],int(entry[1])-1
	# define the left most alignment
	if int(entry[1])-1 < aln.reference_start: left_cigar,left_left,right_cigar,right_left = entry[3],int(entry[1])-1,aln.cigarstring,aln.referece_start
	left_right = Cigar.Cigar(left_cigar).right_position(left_left)-1
	right_right = Cigar.Cigar(right_cigar).right_position(right_left)-1
	return left_left,left_right,right_left,right_right

if __name__=='__main__':
	args = Args()
	err_fh(args.bam)
	err_fh(args.bed)
	bam = pysam.AlignmentFile(args.bam)
	with open(args.bed,'r') as f:
		for l in f:
			r = l.rstrip().split('\t')
			reads={}
			chrom,start,end = r[0],int(r[1]),int(r[2])
			left_region = '{}:{}-{}'.format(chrom,start-args.flank,start+args.flank)
			right_region = '{}:{}-{}'.format(chrom,end-args.flank,end+args.flank)
			reads=fetch_reads(bam,[left_region,right_region]):
			deletion_breakpoint(reads)
			
