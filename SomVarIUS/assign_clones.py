from mixture_model import mixture_model
import numpy as np
import scipy.stats as stats
import argparse, sys


def classify_clones(allele_frequencies, coverage=None, tumor_purity=1.0, gmm=False):
    
    clone_labels = np.empty(len(allele_frequencies), dtype='S9')
    
    if gmm == True:
        best_model, labels, bic = mixture_model(allele_frequencies, 10, p_mean=tumor_purity/2, p_std=0.1, quiet=True)
        c = best_model.components
        clonal_label = 0
        clonal_mean = 0
        for i in xrange(len(c)):
            mean = c[i][0].mu
            #std = c[i][0].sigma
            if mean > clonal_mean:
                clonal_mean = mean
                clonal_label = i
                
        clone_labels[labels == clonal_label] = 'clone'
        clone_labels[labels != clonal_label] = 'sub-clone'
        
        return clone_labels
    
    else:
        #binomial = np.vectorize(stats.binom_test)
        binomial = np.vectorize(stats.binom.cdf)
        alt_alleles = np.array(allele_frequencies*coverage, dtype=int)
        pvals = binomial(alt_alleles, coverage, tumor_purity/2)
        clone_labels[pvals > 0.05] = 'clone'
        clone_labels[pvals <= 0.05] = 'sub-clone'
        
        return clone_labels

def clones(args):
    vcf = args.vcf
    tumor_purity = args.t
    gmm = args.gmm
    out = sys.stdout
    fmt = '{chrom}\t{pos}\t{status}\n'
    
    allele_frequencies = []
    coverage = []
    chroms = []
    pos = []
    for line in open(vcf, 'r'):
        if line[0] != '#':
            fields = line.strip().split('\t')
            info = fields[7].split(';')
            allele_frequencies.append(float(info[1].split('=')[1]))
            coverage.append(float(info[5].split('=')[1]))
            chrom.append(fields[0])
            pos.append(fields[1])
            
    allele_frequencies = np.array(allele_frequencies)
    coverage = np.array(coverage)
    
    clone_labels = classify_clones(allele_frequencies, coverage, tumor_purity)
    
    for i in xrange(len(clone_labels)):
        chrom = chroms[i]
        pos = pos[i]
        status = clone_labels[i]
        out.write(fmt.format(**locals()))
    
    out.close()
    

if __name__ == '__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('--vcf', help='vcf file', required=True)
    parser.add_argument('--t', help='tumor purity (default=1.0)', default=1.0, type=float)
    parser.add_argument('--gmm', help='flag to classify by gaussian mixture model (default=False)', default=False, action='store_true')
    
    args=parser.parse_args()
    
    clones(args)