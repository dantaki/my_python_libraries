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
		parser.add_argument('-o', help='output file',required=False,type=str,default="breaks.out")
		args = parser.parse_args()
		self.bam,self.bed,self.out = args.i,args.b,args.o
		self.flank = args.n
def err_fh(fh):
	if not os.path.isfile(fh): 
		sys.stderr.write('FATAL ERROR: {} NOT FOUND\n'.format(fh))
		sys.exit(1)

def fetch_reads(bam,region):
	for aln in bam.fetch(region=region):
		left_aln,right_aln = None,None #left and right alignments for the breakpoint 
		if not aln.has_tag('SA'):continue
		prm_mapq = aln.mapping_quality
		parse_supplementary_alignments(aln,'DEL')
def parse_supplementary_alignments(aln,svtype):
	"""
	check each supplementary alignment for SV evidence
	"""
	prm_strand = '+' # forward
	if aln.is_reverse: prm_strand='-' #reverse
	sec_alignments = str(aln.get_tag('SA')).split(';')
	if not sec_alignments[-1].startswith('SA'): del sec_alignments[-1] # removes empty entry 
	print len(sec_alignments)
	for sec_aln in sec_alignments:
		entry = sec_aln.split(',')
		sec_contig, sec_left_pos, sec_strand, sec_cigar, sec_mapq = entry[0:5]
		# skip if alignments are on not the same strand
		if prm_strand != sec_strand: continue
		# skip if alignments are not on the same chromosome 
		if sec_contig != aln.reference_name:continue
		if svtype=='DEL': is_deletion(aln,entry)
def is_deletion(aln,entry):
	"""
	check if the two alignments support a deletion

	input: primary alignment object from AlignmentFile().fetch(), list of a SA tag entry
	"""
	left_cigar,left_left,right_cigar,right_left = aln.cigarstring,aln.reference_start,entry[3],int(entry[1])-1
	# define the left most alignment
	if int(entry[1])-1 < aln.reference_start: left_cigar,left_left,right_cigar,right_left = entry[3],int(entry[1])-1,aln.cigarstring,aln.referece_start
	left_right = Cigar.Cigar(left_cigar).right_position(left_left)
	right_right = Cigar.Cigar(right_cigar).right_position(right_left)
	print left_left,left_right,'---',right_left,right_right

if __name__=='__main__':
	args = Args()
	err_fh(args.bam)
	err_fh(args.bed)
	bam = pysam.AlignmentFile(args.bam)
	with open(args.bed,'r') as f:
		for l in f:
			r = l.rstrip().split('\t')
			chrom,start,end = r[0],int(r[1]),int(r[2])
			left_region = '{}:{}-{}'.format(chrom,start-args.flank,start+args.flank)
			right_region = '{}:{}-{}'.format(chrom,end-args.flank,end+args.flank)
			fetch_reads(bam,left_region)