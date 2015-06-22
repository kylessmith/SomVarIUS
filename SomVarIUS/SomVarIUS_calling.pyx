# cython: profile=False
# cython: embedsignature=True


__version__ = '0.15'
#try cimport BedTool, Interval
from collections import OrderedDict
import scipy.stats
from mix_bins import find_bins
import cPickle
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport isnan
from cpython cimport PyBytes_FromStringAndSize
from pysam.chtslib cimport bam1_t, bam_pileup1_t
from pysam.cfaidx cimport Fastafile
from pysam.calignmentfile cimport PileupColumn
from pysam.csamfile cimport Samfile, \
    pysam_bam_get_seq, pysam_bam_get_qual

from scipy.special import gammaln
from scipy.misc import factorial
import numpy as np
cimport numpy as np


cdef tuple fit_beta_binomial(np.ndarray[np.double_t, ndim=1] xi,
 	 						np.ndarray[np.double_t, ndim=1] ni):
							
	cdef np.ndarray[np.float64_t, ndim=1] pi, wi
	cdef int i
	cdef long double w, p, q, s, mu, gamma, old_gamma
	cdef long double rho, beta, alpha

	pi=xi/ni

	# initializing weights and gamma. The weights will be updated recursively later on.
	wi = np.repeat(1.0, len(xi))
	gamma=1

	for i in xrange(1,100):
		w = sum(wi)
		p = sum(wi*pi)/w
		q = (1-p)

		S = sum(wi*(pi-p)**2)

		mu = p
		old_gamma = gamma
		gamma = (S-p*q*(sum((wi/ni)*(1-wi/w))))/(p*q*(sum(wi*(1-wi/w))-sum((wi/ni)*(1-wi/w))))
		wi=ni/(1+gamma*(ni-1))
		if(abs(old_gamma-gamma)<0.0001):
			break

	rho=theta=gamma/(1-gamma)
	beta=(1-mu)*(1/gamma-1)
	alpha=mu*beta/(1-mu)

	return alpha, beta


cdef extern from "math.h":
	long double pow(long double x, long double y)
	long double exp(long double x)
	long double fmax(long double x, long double y)
	long double sqrt(long double x)

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


cdef long double log_dcm(np.ndarray[np.int_t, ndim=1] counts,
 					np.ndarray[np.double_t, ndim=1] alpha):
	N = sum(counts)
	A = sum(alpha)
	return gammaln(N+1) - sum(gammaln(counts+1)) + gammaln(A) - gammaln(N + A) + sum(gammaln(counts+alpha) - gammaln(alpha))


cdef long double log_betabin(int k, int N,
 							long double a, long double b):
	return log_dcm(np.array([k,N-k]),np.array([a,b]))


cdef long double betabin_cdf(int k, int N,
 							long double a, long double b):
	cdef int x

	if k > 0.5 * N:
		p = 1. - sum([np.exp(log_betabin(x,N,a,b)) for x in xrange(k+1,N)])
	else:
		p = sum([np.exp(log_betabin(x,N,a,b)) for x in xrange(k+1)])
	return p


cdef list generate_chrom_dist(chroms, lengths):

	cdef list chrom_regions = []
	cdef str chrom
	cdef int start, end
	for i in xrange(len(chroms)):
		chrom = chroms[i]
		start = 1
		end = lengths[i]
		chrom_regions.append([chrom, start, end])
	return chrom_regions


cdef tuple track_copy_pos(str chrom, int pos, list copy_list, set seen_chroms, int copy_list_pos):

	cdef bint in_copy = False
	cdef str copy_chrom
	cdef int copy_start, copy_end

	for i in xrange(copy_list_pos, len(copy_list)):
		copy_chrom = copy_list[i][0]
		copy_start = int(copy_list[i][1])
		copy_end = int(copy_list[i][2])
		if copy_chrom == chrom:
			seen_chroms.add(chrom)
			if pos >= copy_start and pos <= copy_end:
				in_copy = True
				return copy_chrom, copy_start, copy_end, i, seen_chroms, in_copy
			elif pos < copy_start:
				return copy_chrom, copy_start, copy_end, i, seen_chroms, in_copy
			else:
				pass
		elif chrom in seen_chroms:
			return copy_chrom, copy_start, copy_end, i, seen_chroms, in_copy
		else:
			pass

	copy_chrom = copy_list[copy_list_pos][0]
	copy_start = int(copy_list[copy_list_pos][1])
	copy_end = int(copy_list[copy_list_pos][2])
	in_copy = False
	return copy_chrom, copy_start, copy_end, copy_list_pos, seen_chroms, in_copy


def get_dist(copy_record, Samfile samfile, dict ba_to_bed,
 			int min_reads, int min_support, int min_mapq,
			int min_baseq, dict genome, bint copy_bed=False):

	cdef _CountAllele A, T, G, C
	cdef str chrom, alnbase, maj_allele, alt_allele, base
	cdef int pos, i, maj_count, alt_count
	cdef PileupColumn col
	cdef bam_pileup1_t ** plp
	cdef bam_pileup1_t * read
	cdef bam1_t * aln
	cdef tuple alleles
	cdef bint is_reverse = False
	cdef bint is_proper_pair
	cdef int mapq
	cdef double baseq = 0
	cdef tuple counts, labels, copy_region
	cdef list copy_list = []
	cdef list binned_data
	cdef int n
	cdef long double alpha, beta, mean, total_counts
	cdef str dbsnp_maj_allele
	cdef frozenset dbsnp_alt_alleles
	cdef tuple previous_copy_region = tuple(['None', 'None', 'None'])

	#record = defaultdict(lambda : defaultdict(double))
	record = OrderedDict()
	it = samfile.pileup(reference=None, start=None, end=None,
	truncate=False, max_depth=8000)

	for col in it:
		chrom = samfile.getrname(col.tid)
		pos = col.pos + 1
		_init_count(&A)
		_init_count(&T)
		_init_count(&G)
		_init_count(&C)

		try:
			alleles = copy_record[(ba_to_bed[chrom],pos)][3]
			dbsnp_maj_allele = alleles[0]
			dbsnp_alt_alleles = alleles[1]
			n = col.n
			plp = col.plp
			for i from 0<=i < n:
				read = &(plp[0][i])
				if read.indel == 0 and read.is_del == False:
					aln = read.b
					alnbase = _get_seq_base(aln, read.qpos)
					flag = aln.core.flag
					is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
					mapq = aln.core.qual
					baseq = pysam_bam_get_qual(aln)[read.qpos]
					if mapq >= min_mapq and baseq >= min_baseq and is_proper_pair:
						if alnbase == b'A':
							_incr_count(&A, is_reverse, baseq, mapq)
						elif alnbase == b'T':
							_incr_count(&T, is_reverse, baseq, mapq)
						elif alnbase == b'G':
							_incr_count(&G, is_reverse, baseq, mapq)
						elif alnbase == b'C':
							_incr_count(&C, is_reverse, baseq, mapq)
			if A.pres+T.pres+G.pres+C.pres > 1:
				counts, labels = get_highest(A.count,T.count,G.count,C.count)
				obs_alleles = set([labels[i] for i in xrange(len(labels)) if counts[i]>0])
				#print
				#print "obs", obs_alleles
				dbsnp_obs_overlap = obs_alleles & dbsnp_alt_alleles
				#print "dbsnp", dbsnp_obs_overlap
				#print "dbsnp", dbsnp_maj_allele
				if len(dbsnp_obs_overlap) > 0 and dbsnp_maj_allele in obs_alleles:
					copy_region = copy_record[(ba_to_bed[chrom],pos)][0:3]
					try:
						record[copy_region]
					except KeyError:
						record[copy_region]={'total_alleles':[], 'alt_alleles':[], 'pos':[]}
						if previous_copy_region[1] != 'None':
							if len(record[previous_copy_region]['alt_alleles']) >= 4:
								chrom, start, end = previous_copy_region
								if copy_bed:
									alpha, beta = fit_beta_binomial(np.array(record[previous_copy_region]['alt_alleles'],dtype=np.float64),
																									np.array(record[previous_copy_region]['total_alleles'],dtype=np.float64))
									mean = alpha / (alpha+beta)
									copy_list.append((chrom, start, end, alpha, beta, mean))
								else:
									binned_data = find_bins(record[previous_copy_region]['alt_alleles'],
												record[previous_copy_region]['total_alleles'],
												record[previous_copy_region]['pos'], chrom, genome[chrom], max_components=20)
								
									copy_list += binned_data
							del record[previous_copy_region]
						previous_copy_region = copy_region
					#print chrom, pos
					#print copy_region
					total_counts = float(A.count+T.count+C.count+G.count)
					if total_counts >= min_reads:
						for base in dbsnp_obs_overlap:
							if base == 'A' and A.count >= min_support and total_counts-A.count >= min_support:
								record[copy_region]['alt_alleles'].append(A.count)
								record[copy_region]['total_alleles'].append(total_counts)
								record[copy_region]['pos'].append(pos)
								break
							elif base == 'T' and T.count >= min_support and total_counts-T.count >= min_support:
								record[copy_region]['alt_alleles'].append(T.count)
								record[copy_region]['total_alleles'].append(total_counts)
								record[copy_region]['pos'].append(pos)
								break
							elif base == 'G' and G.count >= min_support and total_counts-G.count >= min_support:
								record[copy_region]['alt_alleles'].append(G.count)
								record[copy_region]['total_alleles'].append(total_counts)
								record[copy_region]['pos'].append(pos)
								break
							elif base == 'C' and C.count >= min_support and total_counts-C.count >= min_support:
								record[copy_region]['alt_alleles'].append(C.count)
								record[copy_region]['total_alleles'].append(total_counts)
								record[copy_region]['pos'].append(pos)
								break

		except KeyError:
			pass
	
	try:
		if len(record[previous_copy_region]['alt_alleles']) >= 4:
			if copy_bed:
				alpha, beta = fit_beta_binomial(np.array(record[previous_copy_region]['alt_alleles'],dtype=np.float64),
																				np.array(record[previous_copy_region]['total_alleles'],dtype=np.float64))
				mean = alpha / (alpha+beta)
				copy_list.append((chrom, start, end, alpha, beta, mean))
			else:
				binned_data = find_bins(record[previous_copy_region]['alt_alleles'],
							record[previous_copy_region]['total_alleles'],
							record[previous_copy_region]['pos'], chrom, genome[chrom], max_components=20)
			
				copy_list += binned_data
		del record[previous_copy_region]
	except KeyError:
		pass

	return copy_list


cdef tuple filter_regions(str dbsnp_bed_name, list regions):

	cdef set seen_chroms = set()
	cdef tuple alleles
	cdef int copy_list_pos = 0
	cdef str line, dbsnp, chrom
	cdef int pos

	copy_record = {}

	for dbsnp in open(dbsnp_bed_name,'r'):
		fields = dbsnp.strip().split('\t')
		chrom = fields[0]
		pos = int(fields[1])
		copy_chrom, copy_start, copy_end, copy_list_pos, seen_chroms, in_copy = track_copy_pos(chrom, pos,
																				regions, seen_chroms, copy_list_pos)
		if in_copy:
			alleles = (fields[3], frozenset(fields[4].split(',')))
			copy_record[(chrom,pos)] = (copy_chrom,copy_start,copy_end,alleles)
		else:
			continue

	return copy_record, seen_chroms


cdef long double calc_se_pval(list copy_list, int pos, str chrom, int alt_count, int total_alleles,
 							dict ba_to_bed, int copy_list_pos, bint in_copy, dict chrom_dist,
 							bint binom):

	cdef list copy_hits
	cdef long double alpha, beta, mean

	if binom:
		pvalue = scipy.stats.binom.cdf(alt_count, total_alleles, 0.5)

	elif in_copy:
		fields = copy_list[copy_list_pos]
		alpha = fields[3]
		beta = fields[4]
		if alpha < 1e15 and beta < 1e15:
			pvalue = betabin_cdf(alt_count, total_alleles, alpha, beta)
		else:
			mean = fields[5]
			pvalue = scipy.stats.binom.cdf(alt_count, total_alleles, mean)
	else:
		try:
			alpha = chrom_dist[chrom][0]
			beta = chrom_dist[chrom][1]
			if alpha < 1e15 and beta < 1e15:
				pvalue = betabin_cdf(alt_count, total_alleles, alpha, beta)
			else:
				mean = chrom_dist[chrom][2]
				pvalue = scipy.stats.binom.cdf(alt_count, total_alleles, mean)
		except:
			pvalue = scipy.stats.binom.cdf(alt_count, total_alleles, 0.5)

	return 1-pvalue


cdef inline tuple chernoff(double m, double d, int alt_count, int total_alleles,
 					double min_se, double min_hetero, int pos, str chrom,
 					bint copy_filter, list copy_list, dict ba_to_bed,
 					bint in_copy, int copy_list_pos, dict chrom_dist,
 					bint binom):

	cdef long double prob_se, prob_not_se, prob_not_hetero, pvalue

	prob_se = pow((exp(d) / pow((1+d), (1+d))), m)
	prob_se=fmax(1e-15, prob_se)
	prob_not_hetero = calc_se_pval(copy_list, pos, chrom, alt_count, total_alleles,
										ba_to_bed, copy_list_pos, in_copy, chrom_dist,
										binom)
	prob_not_se = 1.0 - prob_se
	pvalue = 1.0-(prob_not_se * prob_not_hetero)
	if prob_not_se <= min_se or prob_not_hetero <= min_hetero:
		pvalue = 1.0
	return pvalue, prob_not_se, prob_not_hetero


cdef inline tuple multi_alleles(str alt_allele, str third_allele, str fourth_allele,
 						_CountAllele A, _CountAllele T,
 						_CountAllele G, _CountAllele C,
 						int alt_count, int third_count,
						int fourth_count):

	cdef double alt_quality, third_quality, fourth_quality

	if alt_allele == 'A':
		alt_quality = A.qual / A.count
	elif alt_allele == 'T':
		alt_quality = T.qual / T.count
	elif alt_allele == 'G':
		alt_quality = G.qual / G.count
	elif alt_allele == 'C':
		alt_quality = C.qual / C.count

	if third_allele == 'A':
		third_quality = A.qual / A.count
	elif third_allele == 'T':
		third_quality = T.qual / T.count
	elif third_allele == 'G':
		third_quality = G.qual / G.count
	elif third_allele == 'C':
		third_quality = C.qual / C.count

	if third_count == fourth_count:
		if fourth_allele == 'A':
			fourth_quality = A.qual / A.count
		elif fourth_allele == 'T':
			fourth_quality = T.qual / T.count
		elif fourth_allele == 'G':
			fourth_quality = G.qual / G.count
		elif fourth_allele == 'C':
			fourth_quality = C.qual / C.count

		if alt_quality >= fourth_quality and alt_quality >= third_quality:
			return alt_allele, alt_count
		elif third_quality >= alt_quality and third_quality >= fourth_quality:
			return third_allele, third_count
		elif fourth_quality >= third_quality and fourth_quality >= alt_quality:
			return fourth_allele, fourth_count

	else:
		if third_quality > alt_quality:
			return third_allele, third_count
		else:
			return alt_allele, alt_count


cdef inline double rna_seq_filter(str chrom, int pos, str alt_allele, Samfile rna_seq_file):
	cdef double rna_seq_support = float("NaN")
	cdef bam_pileup1_t ** rna_plp
	cdef bam_pileup1_t * rna_read
	cdef bam1_t * rna_aln
	cdef PileupColumn col
	cdef int rna_pos
	cdef str rna_alnbase
	cdef int i
	cdef int n

	it = rna_seq_file.pileup(reference=chrom, start=pos, end=pos+1,
						truncate=False, max_depth=8000)
	for col in it:
		# initialise variables
		n = col.n
		rna_pos = col.pos
		rna_plp = col.plp
		if rna_pos == pos:
			for i from 0<=i < n:
				rna_read = &(rna_plp[0][i])
				rna_aln = rna_read.b
				rna_alnbase = _get_seq_base(rna_aln, rna_read.qpos)
				if rna_alnbase == alt_allele:
					rna_seq_support =  1
				else:
					rna_seq_support =  0
			return rna_seq_support
	return rna_seq_support


cdef int hapmap_filter(str chrom, int pos, hapmap_record):
	cdef int hapmap_present = 0
	hapmap_present = 0
	try:
		hapmap_present = hapmap_record[chrom][pos]
	except:
		pass
	return hapmap_present


cdef inline tuple get_highest(int A,int T,int G,int C):

	cdef str lA, lT, lG, lC

	lA, lT, lG, lC = "A", "T", "G", "C"
	if A < T:
		A,T,lA,lT = T,A,lT,lA
	if G < C:
		G,C,lG,lC = C,G,lC,lG
	if A < G:
		A,G,lA,lG = G,A,lG,lA
	if T < C:
		T,C,lT,lC = C,T,lC,lT
	if T < G:
		T,G,lT,lG = G,T,lG,lT
	return (A,T,G,C), (lA,lT,lG,lC)


cdef dict samfile_ref_dicts(samfile_refs, other_refs):
	cdef str ref
	cdef dict samfile_to_other_ref={}
	for ref in samfile_refs:
		if "chr"+ref in other_refs:
			samfile_to_other_ref[ref]="chr"+ref
		elif ref[3:] in other_refs:
			samfile_to_other_ref[ref]=ref[3:]
		else:
			samfile_to_other_ref[ref]=ref
	return samfile_to_other_ref


def write_mutations(str out_file, str samfile_name, str fafile_name,
						  str rna_seq_file_name='', str hapmap_file='',
                          int min_reads=10, int min_support=4, double min_af=0.05,
						  double min_pvalue=0.0001, double min_f_r=0.05,
						  bint rna_filter=False, bint hap_filter=False,
						  int min_qual=25, double min_se=0.999,
						  double min_hetero=0.95, bint ref_filter=False,
						  uint64_t min_mapq=55, str dbsnp_bed_name='',
						  str copy_bed_name='', bint binom=False,
						  region_chrom=None, region_start=None, region_end=None,
						  int min_baseq=13, str dist='', bint whole=False):

	#out = open(out_file, mode='w')
	fmt = '{chrom}\t{pos}\t.\t{refbase}\t{alt_allele}\t.\tPASS\t'
	fmt+='NS=1;AF={alt_allele_freq};P={pvalue};SP={prob_not_se};'
	fmt+='GP={prob_not_hetero};DP={total_alleles};MQ={mean_mapq};BQ={mean_qual}\t'
	fmt+='GT:GQ:DP\t1|0:60:{total_alleles}\n'
	
	#fmt='{chrom}\t{pos}\t{refbase}\t{maj_allele}\t{alt_allele}'
	#fmt+='\t{alt_allele_freq}\t{pvalue}\t{prob_not_se}\t{prob_not_hetero}'
	#fmt+='\t{total_alleles}\t{alt_count}\t{rvs_ratio}\t{fwd_ratio}\t{mean_qual}\t{mean_mapq}\n'
	cdef Samfile rnq_seq_file
	cdef Samfile samfile
	cdef Fastafile fafile
	cdef dict ba_to_fa, ba_to_rna, ba_to_hapmap, ba_to_bed = {}
	cdef set bed_chroms, seen_chroms = set()
	cdef list copy_list
	cdef int copy_list_pos = 0
	cdef bint in_copy = False
	cdef str copy_chrom
	cdef int copy_start, copy_end

	samfile = Samfile(samfile_name)
	fafile = Fastafile(fafile_name)
	
	ba_to_fa = samfile_ref_dicts(samfile.references, fafile.references)
	
	references = samfile.references
	lengths = samfile.lengths
	genome = {references[ref]:lengths[ref] for ref in xrange(len(references))}
	# if copy number file provided
	if copy_bed_name != '':
		copy_filter = True
		print 'storing dbsnp/copy number intersection in memory...'
		copy_regions = []
		for line in open(copy_bed_name,'r'):
			copy_regions.append(line.strip().split('\t'))
		copy_record, bed_chroms = filter_regions(dbsnp_bed_name, copy_regions)
		ba_to_bed = samfile_ref_dicts(samfile.references, bed_chroms)
		print 'calculating allele frequency distribution...'
		copy_list = get_dist(copy_record, samfile, ba_to_bed, min_reads, min_support, min_mapq, min_baseq, genome, copy_bed=True)
		if dist != '':
			dist_out = open(dist,"w")
			for line in copy_list:
				dist_out.write('\t'.join(map(str,line))+'\n')
		print 'done calculating copy number distributions'
	
		del copy_record
		print 'storing chromosome/dbsnp intersection in memory...'
		chrom_regions = generate_chrom_dist(samfile.references, samfile.lengths)
		chrom_record, all_bed_chroms = filter_regions(dbsnp_bed_name, chrom_regions)
		print 'calculating chromosome allele frequency distribution...'
		chrom_list = get_dist(chrom_record, samfile, ba_to_bed, min_reads, min_support, min_mapq, min_baseq, genome, copy_bed=True)
		if dist != '':
			for line in chrom_list:
				dist_out.write('\t'.join(map(str,line))+'\n')
			dist_out.close()
		chrom_dist = {field[0]:(field[3],field[4], field[5]) for field in chrom_list}
		del chrom_record, all_bed_chroms, chrom_list
	
		copy_chrom = copy_list[0][0]
		copy_start = int(copy_list[0][1])
		copy_end = int(copy_list[0][2])
		print 'done calculating chromosome distributions'
		print
	
	# if no copy number file provided
	else:
		copy_filter = True
		print 'storing chromosome/dbsnp intersection in memory...'
		chrom_regions = generate_chrom_dist(samfile.references, samfile.lengths)
		chrom_record, bed_chroms = filter_regions(dbsnp_bed_name, chrom_regions)
		ba_to_bed = samfile_ref_dicts(samfile.references, bed_chroms)
		print 'calculating chromosome allele frequency distribution...'
		if whole:
			copy_list = get_dist(chrom_record, samfile, ba_to_bed, min_reads, min_support, min_mapq, min_baseq, genome, copy_bed=True)
		else:
			copy_list = get_dist(chrom_record, samfile, ba_to_bed, min_reads, min_support, min_mapq, min_baseq, genome, copy_bed=False)
		if dist != '':
			dist_out = open(dist,"w")
			for line in copy_list:
				dist_out.write('\t'.join(map(str,line))+'\n')
			dist_out.close()
		chrom_dist = {field[0]:(field[3],field[4], field[5]) for field in copy_list}
		del chrom_record, chrom_regions
	
		copy_chrom = ''
		copy_start = 0
		copy_end = 0		
	
		print 'done calculating chromosome distributions'
		print
	
	if rna_filter:
		fmt = '{chrom}\t{pos}\t{refbase}\t{maj_allele}\t{alt_allele}'
		fmt+='\t{alt_allele_freq}\t{pvalue}\t{prob_not_se}\t{prob_not_hetero}'
		fmt+='\t{n}\t{alt_count}\t{rvs_ratio}\t{fwd_ratio}\t{mean_qual}\t{mean_mapq}\t{rna_seq}\n'
		rna_seq_file = Samfile(rna_seq_file_name)
		ba_to_rna = samfile_ref_dicts(samfile.references, rna_seq_file.references)
			
	if hap_filter:
		print 'loading germline record...'
		hapmap_record = cPickle.load(open(hapmap_file,'rb'))
		ba_to_hapmap = samfile_ref_dicts(samfile.references, hapmap_record.keys())
		print 'done loading germline record'
		print

	#out.write(fmt.replace('}','').replace('{',''))
	samfile = Samfile(samfile_name)
	fafile = Fastafile(fafile_name)

	it = samfile.pileup(reference=region_chrom, start=region_start, end=region_end,
						truncate=False, max_depth=80000)
			  
	# statically typed variables
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
	cdef int hapmap_present
	cdef double mean_qual
	cdef _CountAllele A, T, G, C
	cdef str maj_allele, alt_allele, third_allele, fourth_allele
	cdef int maj_count, alt_count, third_count, fourth_count
	cdef double fwd_ratio, rvs_ratio
	cdef int total_alleles
	cdef long double pvalue, prob_not_se, prob_not_hetero
	cdef double delta, reverse, forward
	cdef double rna_seq

	print 'calling mutations'

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
		
			# check if hapmap filter is on
			if hap_filter:
				hapmap_present=hapmap_filter(ba_to_hapmap[chrom], pos, hapmap_record)
				# if position is in hapmap file, then skip
				if hapmap_present == 1:
					continue
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
		
			if A.pres+T.pres+G.pres+C.pres > 1:
				counts, labels = get_highest(A.count,T.count,G.count,C.count)
				maj_allele = labels[0]
				maj_count = counts[0]
				alt_allele = labels[1]
				alt_count = counts[1]
				if maj_count != alt_count:
					if A.pres+T.pres+G.pres+C.pres > 2:
						third_allele = labels[2]
						third_count = counts[2]
						fourth_allele = labels[3]
						fourth_count = counts[3]
						if third_count == alt_count:
							alt_allele, alt_count = multi_alleles(alt_allele, third_allele, fourth_allele,
													A, T, G, C, alt_count, third_count, fourth_count)
					if alt_allele == 'A':
						mean_mapq = A.mapqual / float(A.count)
						mean_qual = A.qual / A.count
						fwd_ratio = A.fwd / float(A.fwd+A.rev)
						rvs_ratio = A.rev / float(A.fwd+A.rev)
					elif alt_allele == 'T':
						mean_mapq = T.mapqual / float(T.count)
						mean_qual = T.qual / T.count
						fwd_ratio = T.fwd / float(T.fwd+T.rev)
						rvs_ratio = T.rev / float(T.fwd+T.rev)
					elif alt_allele == 'G':
						mean_mapq = G.mapqual / float(G.count)
						mean_qual = G.qual / G.count
						fwd_ratio = G.fwd / float(G.fwd+G.rev)
						rvs_ratio = G.rev / float(G.fwd+G.rev)
					elif alt_allele == 'C':
						mean_mapq = C.mapqual / float(C.count)
						mean_qual = C.qual / C.count
						fwd_ratio = C.fwd / float(C.fwd+C.rev)
						rvs_ratio = C.rev / float(C.fwd+C.rev)

					if alt_count >= min_support and fwd_ratio >= min_f_r and rvs_ratio >= min_f_r and mean_qual >= min_qual and mean_mapq >= min_mapq:
						total_baseq = fmax(1e-25,total_baseq)
						delta = (alt_count/total_baseq) - 1
						total_alleles = A.count+T.count+C.count+G.count
						if copy_filter:
							if copy_chrom == chrom and pos >= copy_start and pos <= copy_end:
								in_copy = True
							elif copy_list_pos >= (len(copy_list)-1):
								in_copy = False
							else:
								copy_chrom, copy_start, copy_end, copy_list_pos, seen_chroms, in_copy = track_copy_pos(ba_to_bed[chrom],
																										pos, copy_list, seen_chroms,
																										copy_list_pos)
						else:
							in_copy = False
						pvalue, prob_not_se, prob_not_hetero = chernoff(total_baseq, delta, alt_count,
																total_alleles, min_se, min_hetero, pos,
																chrom, copy_filter, copy_list, ba_to_bed,
																in_copy, copy_list_pos, chrom_dist, binom)
						if pvalue <= min_pvalue:
							alt_allele_freq = alt_count / float(total_alleles)
							if alt_allele_freq >= min_af:
								# reference base
								prob_not_se = 1 - prob_not_se
								prob_not_hetero = 1 - prob_not_hetero
								refbase = fafile\
								    .fetch(reference=ba_to_fa[chrom], start=col.pos, end=col.pos+1)\
								    .upper()
								if rna_filter:
									rna_seq = rna_seq_filter(ba_to_rna[chrom],col.pos,alt_allele,rna_seq_file)
								if ref_filter:
									if refbase != alt_allele:
										yield fmt.format(**locals())
								else:
									yield fmt.format(**locals())
	#out.close()


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