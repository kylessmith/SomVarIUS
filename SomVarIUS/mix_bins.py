import numpy as np
import mixture_model
from bb_fit import bb_model
#from pymix import mixture
import argparse
from collections import OrderedDict

def find_bins(alt_freqs, total_freqs, pos, chrom, max_length, max_components=20):
    
    alt_freqs = np.array(alt_freqs)
    total_freqs = np.array(total_freqs, dtype=float)
    pos = np.array(pos)
    #data = mixture.DataSet()
    #data.fromArray(alt_freqs/total_freqs)
    
    allele_freqs = np.array(alt_freqs / total_freqs)

    best_model, labels, best_model_bic = mixture_model.mixture_model(allele_freqs,
                                max_components, p_mean=0.5, p_std=0.1, quiet=True)
        
    model_lookup = {}
    
    for i in np.unique(labels):
        model = bb_model(alt_freqs[labels==i], total_freqs[labels==i])
        model_lookup[i] = model
    
    pos_starts = pos[np.arange(0,len(alt_freqs),5)]
    pos_ends = pos[np.arange(4,len(alt_freqs),5)]
    
    if len(pos_ends) < len(pos_starts):
        pos_ends = np.append(pos_ends, max_length)
    else:
        pos_ends[-1] = max_length
    
    pos_starts[0] = 0
    diffs = pos_starts[1:] - pos_ends[:-1]
    s_diffs = diffs/2
    e_diffs = diffs-s_diffs
    pos_starts = np.append(0, pos_starts[1:]-s_diffs)
    pos_ends = np.append(e_diffs+pos_ends[:-1], max_length)
    
    bins = np.reshape(labels[:-(len(labels)%5)], (-1,5))
    
    label_bins = []
    for i in xrange(len(bins)):
        l = np.argmax(np.bincount(bins[i]))
        start = pos_starts[i]
        if i == len(bins)-1:
            end = max_length
        else:
            end = pos_ends[i]
        label_bins.append( (chrom,start,end,model_lookup[l].alpha, 
                            model_lookup[l].beta, model_lookup[l].mean) )
        
    return label_bins
    
def read_germ(fname):
    
    chrom_pos = OrderedDict()
    chrom_alt = OrderedDict()
    chrom_total = OrderedDict()
    
    for line in open(fname,'r'):
        fields = line.strip().split('\t')
        chrom = fields[0]
        if chrom not in chrom_pos:
            chrom_pos[chrom] = []
            chrom_alt[chrom] = []
            chrom_total[chrom] = []
        
        chrom_pos[chrom].append(int(fields[1]))
        chrom_alt[chrom].append(int(fields[2]))
        chrom_total[chrom].append(float(fields[3]))
        
    for chrom in chrom_pos:
        chrom_pos[chrom] = np.array(chrom_pos[chrom])
        chrom_alt[chrom] = np.array(chrom_alt[chrom])
        chrom_total[chrom] = np.array(chrom_total[chrom])
        
    return chrom_pos, chrom_alt, chrom_total


def read_genome(genome):
    
    chrom_lengths = {}
    
    for line in open(genome,"r"):
        fields = line.strip().split('\t')
        chrom_lengths[fields[0]] = int(fields[1])
    
    return chrom_lengths


def chrom_bins(chrom_lengths, chrom_alt, chrom_total, chrom_pos, max_components):
    
    binned_chroms = []
    
    for chrom in chrom_pos:
        label_bins = find_bins(chrom_alt[chrom], chrom_total[chrom], chrom_pos[chrom],
                                            chrom, chrom_lengths[chrom], max_components)
                                            
        binned_chroms += label_bins
        
    return binned_chroms


def main():
    
    parser=argparse.ArgumentParser()
    parser.add_argument('--i', help='input germline file', required=True)
    parser.add_argument('--g', help='input genome file', required=True)
    parser.add_argument('--max', help='max number of components', default=20, type=int)
    args=parser.parse_args()
    
    fname = args.i
    genome = args.g
    max_components = args.max
    
    chrom_lengths = read_genome(genome)
    chrom_pos, chrom_alt, chrom_total = read_germ(fname)
    binned_chroms = chrom_bins(chrom_lengths, chrom_alt, chrom_total, chrom_pos, max_components)

    
if __name__ == "__main__":
    main()