#!/usr/bin/python
"""
Generating permutation data
============================
"""

def run():
	"""
	This is a simple Python script for generating permuted mutation data. The
	default parameter values are for permuting the GBM(2008) mutation data in
	the examples/mutation_data/GBM_2008 directory, which is a small dataset.
	The main utility of this script is that it is much faster to generate
	permutation data before running the permutation test.

	On my machine, it takes ~3 seconds to permute the mutation data *once*.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../'))
	import multi_dendrix.permute.mutation_data as P

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	MUTATION_MATRIX="../mutation_data/GBM_2008/GBM_2008.m2"
	CUTOFF=2
	NUM_MATRICES=1
	OUTPUT_DIR="./sample_output/permuted_matrices"
	GENE_BLACKLIST="../mutation_data/fishy_genes.glst"
	GENE_WHITELIST=None
	PATIENT_WHITELIST=None
	PATIENT_BLACKLIST=None

	args = [ "-m", MUTATION_MATRIX, "-c", CUTOFF, "-o", OUTPUT_DIR, "-v" ]

	if GENE_BLACKLIST: args += [ "-bg", GENE_BLACKLIST ]
	if GENE_WHITELIST: args += [ "-g", GENE_WHITELIST ]
	if PATIENT_BLACKLIST: args += [ "-bp", PATIENT_BLACKLIST ]
	if PATIENT_WHITELIST: args += [ "-p", PATIENT_WHITELIST ]
	if NUM_MATRICES: args += [ "-n", NUM_MATRICES ]

	# Use the permute_mutation_data module to generate the permuted matrices
	P.run(P.parse_args(map(str, args)))

if __name__ == "__main__": run()
