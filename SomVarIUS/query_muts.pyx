from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from cpython cimport PyBytes_FromStringAndSize
from pysam.chtslib cimport bam1_t, bam_pileup1_t
from pysam.calignmentfile cimport PileupColumn
from pysam.csamfile cimport Samfile, \
    pysam_bam_get_seq, pysam_bam_get_qual


## These are bits set in the flag.
## have to put these definitions here, in csamtools.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
DEF BAM_FPAIRED       =1
## @abstract the read is mapped in a proper pair */
DEF BAM_FPROPER_PAIR  =2
## @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
DEF BAM_FUNMAP        =4
## @abstract the mate is unmapped */
DEF BAM_FMUNMAP       =8
## @abstract the read is mapped to the reverse strand */
DEF BAM_FREVERSE      =16
## @abstract the mate is mapped to the reverse strand */
DEF BAM_FMREVERSE     =32
## @abstract this is read1 */
DEF BAM_FREAD1        =64
## @abstract this is read2 */
DEF BAM_FREAD2       =128
## @abstract not primary alignment */
DEF BAM_FSECONDARY   =256
## @abstract QC failure */
DEF BAM_FQCFAIL      =512
## @abstract optical or PCR duplicate */
DEF BAM_FDUP        =1024


cdef char* CODE2CIGAR= "MIDNSHP=X"

cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

cdef extern from "math.h":
	long double fmax(long double x, long double y)
	long double pow(long double x, long double y)
	long double exp(long double x)
	

def chernoff(double m, double d):

	prob_se = pow((exp(d) / pow((1+d), (1+d))), m)
	prob_se=fmax(1e-15, prob_se)

	return prob_se


cdef dict read_mutation_bed(str mut_fn):
	
	cdef dict possible_muts = {}
	cdef list fields
	cdef str line
	
	for line in open(mut_fn,'r'):
		fields = line.strip().split('\t')
		possible_muts[(fields[0], int(fields[1]))] = (fields[3], fields[4])
		
		
	return possible_muts


def query_mutations(muts_fn, samfile_name, min_mapq, min_baseq, min_reads, min_support,
					out_fn):
	
	
	#statically type variables
	cdef Samfile samfile
	cdef int current_pos = 0
	cdef set seen_chroms = set()
	cdef dict possible_muts
	cdef bint in_germline
	cdef bam_pileup1_t ** plp
	cdef bam_pileup1_t * read
	cdef bam1_t * aln
	cdef uint32_t flag
	cdef bint is_reverse, is_proper_pair
	cdef int i  # loop index
	cdef int n  # total number of reads in column
	cdef double baseq
	cdef double total_baseq = 0.0
	cdef PileupColumn col
	cdef double alt_allele_freq
	cdef str chrom
	cdef str alnbase
	cdef int pos
	cdef double mean_qual
	cdef _CountAllele A, T, G, C
	cdef _CountAllele alt_mut
	cdef str alt_allele
	cdef double fwd_ratio, rvs_ratio
	cdef int total_alleles
	cdef long double prob_se
	cdef double delta, reverse, forward
	
	
	fmt = '{chrom}\t{pos}\t.\t{refbase}\t{alt_allele}\t.\tPASS\t'
	fmt+='NS=1;AF={alt_allele_freq};P={pvalue};SP={prob_se};'
	fmt+='GP={prob_germline};DP={total_reads};MQ={mean_mapq};BQ={mean_qual}\t'
	fmt+='GT:GQ:DP\t1|0:60:{total_reads}\n'
	
	possible_muts = read_mutation_bed(muts_fn)
	
	out = open(out_fn,'a')
	samfile = Samfile(samfile_name)
	it = samfile.pileup(reference=None, start=None, end=None,
							truncate=False, max_depth=80000)
	
	for col in it:
		# initialise variables
		n = col.n
		if n >= min_reads:
			plp = col.plp
			total_baseq = 0.0
			_init_count(&A)
			_init_count(&T)
			_init_count(&G)
			_init_count(&C)

			# get chromosome name and position
			chrom = samfile.getrname(col.tid)
			pos = col.pos + 1
	
			# loop over reads, extract what we need
			for i from 0<=i < n:
				read = &(plp[0][i])
			
				if read.indel == 0 and read.is_del == False:
					aln = read.b
					flag = aln.core.flag
					mapq = aln.core.qual
					is_reverse = <bint>(flag & BAM_FREVERSE)
					is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
					alnbase = _get_seq_base(aln, read.qpos)
					baseq = pysam_bam_get_qual(aln)[read.qpos]
					total_baseq += pow(10.0, (-baseq/10.0))
				
					if mapq >= min_mapq and baseq >= min_baseq and is_proper_pair:
						if alnbase == b'A':
							_incr_count(&A, is_reverse, baseq, mapq)
						elif alnbase == b'T':
							_incr_count(&T, is_reverse, baseq, mapq)
						elif alnbase == b'G':
							_incr_count(&G, is_reverse, baseq, mapq)
						elif alnbase == b'C':
							_incr_count(&C, is_reverse, baseq, mapq)
			#determine if more than 1 allele is observed
			if A.pres+T.pres+G.pres+C.pres > 1:
				total_reads = float(A.count+T.count+C.count+G.count)
				#check if position is in possible mutations
				try:
					possible_muts[(chrom, pos)]
					#check if there is more than one mutation present
					if len(possible_muts[(chrom, pos)][1]) > 1:
						possible_mut = possible_muts[(chrom, pos)][1].split(',')
						alt_mut = eval(possible_mut[0])
						alt_allele = possible_mut[0]
						for i in xrange(1, len(possible_mut)):
							if eval(possible_mut[i])['count'] > alt_mut.count:
								alt_mut = eval(possible_mut[i])
								alt_allele = possible_mut[i]
					else:
						alt_mut = eval(possible_muts[(chrom, pos)][1])
						alt_allele = possible_muts[(chrom, pos)][1]
					
					if eval(possible_muts[(chrom, pos)][0])['count'] > 0 and alt_mut.count > 0:
					
						if alt_mut.count >= min_support:
							mean_mapq = alt_mut.mapqual / float(alt_mut.count)
							mean_qual = alt_mut.qual / alt_mut.count
							fwd_ratio = alt_mut.fwd / float(alt_mut.fwd+alt_mut.rev)
							rvs_ratio = alt_mut.rev / float(alt_mut.fwd+alt_mut.rev)
							total_baseq = fmax(1e-25,total_baseq)
							delta = (alt_mut.count/total_baseq) - 1
							prob_se = chernoff(total_baseq, delta)
							refbase = possible_muts[(chrom,pos)][0]
							alt_allele_freq = alt_mut.count / total_reads
							if alt_allele_freq <= 0.4:
								prob_germline = 'SOMATIC/AMBIGUOUS'
							else:
								prob_germline = 'GERMLINE/AMBIGUOUS'
							pvalue = 'nan'
							
							out.write(fmt.format(**locals()))
				
				except KeyError:
					pass
					
	out.close()
							
							
							
							
cdef inline object _get_seq_base(bam1_t *src, uint32_t k):
	cdef uint8_t * p
	cdef char * s

	if not src.core.l_qseq:
		return None

	seq = PyBytes_FromStringAndSize(NULL, 1)
	s   = <char*>seq
	p   = pysam_bam_get_seq(src)

	s[0] = bam_nt16_rev_table[p[k//2] >> 4 * (1 - k%2) & 0xf]

	return seq

cdef struct _CountAllele:
	int fwd, rev, pres, count, mapqual
	double qual

cdef inline _init_count(_CountAllele* c):
	c.qual = c.fwd = c.rev = c.pres = c.count = c.mapqual = 0

cdef inline _incr_count(_CountAllele* c, bint is_reverse, double baseq, int mapq):
	c.count+=1
	c.pres=1
	c.qual+=baseq
	c.mapqual+=mapq
	if is_reverse:
		c.rev += 1
	else:
		c.fwd += 1