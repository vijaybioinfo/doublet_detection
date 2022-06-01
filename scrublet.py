#/usr/bin/env python

###############################
# Doublet detection: scrublet #
###############################

import optparse
parser=optparse.OptionParser(usage="usage: %prog [options] --input counts")
parser.add_option("-i", "--input",
                  help="Input matrix: unnormalized UMI counts.",
                  metavar="FILE")
parser.add_option("-o", "--outdir", default = "./",
                  help="Output: directory.",
                  metavar="FILE")
parser.add_option("-v", "--verbose", action="store_true", default=True,
                  help="Verbose: Show progress.",
                  metavar="True")

# Getting arguments from command line
(opt, args) = parser.parse_args()

# opt.input = "/home/ciro/large/asthma_biopsy/raw/cellranger/count/Biopsy4_Hu_45P_5P_Gex/outs/filtered_feature_bc_matrix"
# opt.outdir = "/home/ciro/large/asthma_biopsy/results/ab_demux/biop1to7_100th/Biopsy4_Hu_45P_5P"

'''
Given a raw (unnormalized) UMI counts matrix counts_matrix with cells as rows
and genes as columns, calculate a doublet score for each cell:
'''

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc
import time

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

os.chdir(opt.outdir)
print("Working at {}".format(os.getcwd()))
print("Content:"); os.listdir("./")
if os.path.isdir(opt.input):
  print("Reading data:", opt.input)

counts_matrix = sc.read_10x_mtx(opt.input,
                                var_names='gene_symbols',
                                cache=True)
genes = counts_matrix.var.index.values
# counts_matrix = scipy.io.mmread(opt.input + '/matrix.mtx.gz').T.tocsc()
# genes = np.array(scr.load_genes(opt.input + '/genes.tsv.gz', delimiter='\t', column=1))

if True:
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))
    print(genes[:5])

start_time = time.time()
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
print("--- %s seconds ---" % round(time.time() - start_time, 2))

np.savetxt("scrublet_scores.csv", doublet_scores, delimiter=",")
np.savetxt("scrublet_predicted_doublets.csv", predicted_doublets, delimiter=",")

scrub.plot_histogram()
plt.savefig('scrublet_histogram.png', dpi=300, bbox_inches='tight')
