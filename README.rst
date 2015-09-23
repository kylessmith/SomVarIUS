Somatic variant identification from unpaired sample

Data Format
===========

Input files mapped reads must be in SAM or BAM format.

Known germline mutation positions to be excluded must be in a pickled dictionary of dictionaries:

	record[(chromosome, position)] = (reference, alternate)
	
A built-in function is included to create this file from a BED file *(see Invocation)*


Known germline mutations to be used for the prior distribution must be in bed format:

	https://raw.githubusercontent.com/kylessmith/SomVarIUS/master/example/chr20_dbsnp.bed
	
The reference fasta must be in FASTA format.


The copy number segment file must be in BED format (first columns are chrome, start, stop)::

    chr1    11873   14409

NOTE: The BAM and BED files are assumed to be *sorted* in the same way, however a
      built-in function is included to do this *(see Invocation)*

Invocation
==========

Running the following command will the available functions::

    $ SomVarIUS -h

Gives::

	positional arguments:
	  {call_mutations,sort,pickle,clones,query_mutations}
	    call_mutations      flag to call mutations
	    sort                flag to sort bam file and bed file by name
	    pickle              flag to store the pickled germline positions from bed
		clones              flag to classify as clone or sub-clone
		query_mutations     flag to query given mutations in the bam

	optional arguments:
	  -h, --help            show this help message and exit

To call mutations the following command will show available inputs::

	$ SomVarIUS call_mutations -h

Gives::

	  --bam BAM             input bam file
	  --ref REF             reference fasta file
	  --out OUT             output file
	  --rna_seq RNA_SEQ     RNA-seq bam file name
	  --germ_pos GERM_POS   pickled hapmap file
	  --dbsnp_bed DBSNP_BED
	                        dbsnp bed file name
	  --copy_bed COPY_BED   copy number bed file name
	  --min_reads MIN_READS
	                        minimum base coverage (default=10)
	  --min_support MIN_SUPPORT
	                        minimum number of reads supported alternate allele
	                        (default=4)
	  --min_af MIN_AF       minimum allele frequency (default=0.05)
	  --min_pvalue MIN_PVALUE
	                        minimum pvalue (default=0.001)
	  --min_fr MIN_FR       minimum reverse/forward read ratio (default=0.05)
	  --min_qual MIN_QUAL   minimum mean quality for alternate allele (default=25)
	  --min_se MIN_SE       minimum probability not sequencing error
	                        (default=0.999)
	  --min_hetero MIN_HETERO
	                        minimum probability not germline (default=0.95)
	  --ref_filter          flag to filter by reference fasta (default=False)
	  --binom               flag to use binomial test instead of beta-binomial
	                        (default=False)
	  --min_mapq MIN_MAPQ   minimum mapping quality (default=55)
	  --min_baseq MIN_BASEQ
	                        minimum base quality (default=13)
	  --chrom CHROM         Chromosome name to look at
	  --start START         starting position
	  --end END             ending position
	  --dist DIST           write the beta binomial parameters to a file

To sort a bam file and bed file by name the following command will show available inputs::

	$ SomVarIUS sort -h

Gives::

	  --bam BAM             input bam file
	  --bam_out BAM_OUT     name of sorted bam file
	  --dbsnp DBSNP         input dbsnp bed file
	  --dbsnp_out DBSNP_OUT
	                        name of sorted dbsnp file

To create a pickle of bed file positions the following command will show available inputs::

	$ SomVarIUS pickle -h

Gives::

	  --dbsnp DBSNP         input dbsnp bed file
	  --dbsnp_out DBSNP_OUT
	                        name of pickled dbsnp file

To classify mutations as clonal or sub-clonal::

	$ SomVarIUS clones -h

Gives::

	  --vcf VCF   vcf file
	  --t T       tumor purity (default=1.0)
	  --gmm       flag to classify by gaussian mixture model (default=False)

To query a list of mutations in bed format (chrom  start  end  ref  alt)::

	$ SomVarIUS query_mutations -h
	
Gives::

	  --bam BAM             input bam file
	  --out OUT             output file
	  --muts MUTS           mutation bed file
	  --min_reads MIN_READS
	                        minimum base coverage (default=10)
	  --min_support MIN_SUPPORT
	                        minimum number of reads supported alternate allele
	                        (default=4)
	  --min_mapq MIN_MAPQ   minimum mapping quality (default=55)
	  --min_baseq MIN_BASEQ
	                        minimum base quality (default=13)
		
QuickStart
==========

If your files are sorted in the same way and you want to call somatic mutations in all chromosomes.

somatic mutations
-----------------
::

	$ SomVarIUS call_mutations \
		--bam test.bam \
		--ref test.fa \
		--out test_output.txt \
		--germ_pos dbsnp_pos.pickle \
		--dbsnp_bed test_dbsnp.bed \
		--ref_filter

The output will be shown in VCF format.

Example
=======

To run the example files, from the examples directory first run::

	$ SomVarIUS pickle \
		--dbsnp chr20_dbsnp.bed
		--dbsnp_out chr20_dbsnp.pickle
		
Then run::

	$ SomVarIUS call_mutations \
		--bam chr20.bam \
		--ref chr20.fa \
		--out chr20.vcf \
		--germ_pos chr20_dbsnp.pickle \
		--dbsnp_bed chr20_dbsnp.bed \
		--dist dist.txt \
		--min_pvalue 0.05 \
		--ref_filter
		
The first time this is run the program will detect the files have not been index and index them.
The results will be in *chr20.vcf* file and the *dist.txt* will have the estimated parameters 
for the fitted beta-binomial distribution. The arguments used are recorded in the *chr20_args.txt*
file.

Installation
============

pip can be used to install by::

	$ pip install SomVarIUS
	
or download from github and run::

	$ python setup.py install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  

    https://continuum.io/downloads

Dependencies can be installed by::

    pip install -r requirements.txt

The program also depends on Samtools which is available from https://github.com/samtools/samtools