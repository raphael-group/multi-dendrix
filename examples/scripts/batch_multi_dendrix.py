#!/usr/bin/python
"""
Batch Multi-Dendrix
============================
"""

def run():
	"""
	This is a simple Python script for running Multi-Dendrix on different 
	datasets. In this example, I run Multi-Dendrix for k=3, t=3 on each
	of the three mutation datasets: GBM(2008), GBM, and BRCA.

	On my machine, it takes ~20.15 seconds to complete.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../multi_dendrix/lib'))
	import multi_dendrix as Multi

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	K=3
	T=3

	# Set up GBM(2008) args
	MUTATION_MATRIX="../mutation_data/GBM_2008/GBM_2008.m2"
	CUTOFF=2
	OUTPUT_FILE="./sample_output/collections/GBM_2008_t3_k3"
	GENE_BLACKLIST="../mutation_data/fishy_genes.glst"

	gbm08_args = [ "-m", MUTATION_MATRIX, "-c", CUTOFF, "-o", OUTPUT_FILE,
	               "-bg", GENE_BLACKLIST ]

	# Set up GBM args
	MUTATION_MATRIX="../mutation_data/GBM/GBM.m2"
	OUTPUT_FILE="./sample_output/collections/GBM_t3_k3"
	GENE_WHITELIST="../mutation_data/GBM/GBM.glst"

	gbm_args = [ "-m", MUTATION_MATRIX, "-o", OUTPUT_FILE,
	             "-g", GENE_WHITELIST ]

	# Set up BRCA args
	MUTATION_MATRIX="../mutation_data/BRCA/BRCA.m2"
	OUTPUT_FILE="./sample_output/collections/BRCA_t3_k3"
	GENE_WHITELIST="../mutation_data/BRCA/BRCA.glst"

	brca_args = [ "-m", MUTATION_MATRIX, "-o", OUTPUT_FILE,
	              "-g", GENE_WHITELIST ]

	# Run Multi-Dendrix for each of the parameter settings
	for args in [gbm08_args, gbm_args, brca_args]:
		Multi.run(Multi.parse_args(map(str, args + ["-k", K, "-t", T])))

if __name__ == "__main__": run()
