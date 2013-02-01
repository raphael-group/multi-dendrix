#!/usr/bin/python
"""
Weight permutation test
============================
"""

def run():
	"""
	This script performs the matrix permutation test on a collection of gene 
	sets. It gives an example of how one could use various functions from 
	different modules in order to evaluate the collection without obtaining
	the collection from the Multi-Dendrix pipeline. The program first loads
	the GBM(2008) mutation data, permutes it the given number of times, and
	then performs the matrix permutation test.

	The default parameters permute the GBM(2008) data *ten times*.
	On my machine, this script runs in 41.2 seconds.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../'))
	import multi_dendrix.permute.mutation_data as P
	from multi_dendrix.evaluate.matrix import matrix_permutation_test

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	COLLECTION = [["CDKN2B", "RB1", "CDK4"], ["CDKN2A", "TP53", "DTX3"]]
	MUTATION_MATRIX="../mutation_data/GBM_2008/GBM_2008.m2"
	GENE_BLACKLIST="../mutation_data/fishy_genes.glst"
	CUTOFF=2
	ALPHA=1.0
	NUM_PERMUTATIONS=10
	T, K_MIN, K_MAX, OVERLAPS, SETS_PER_GENE  = 2, 3, 3, 0, 1

	# Load the mutation data
	import multi_dendrix as Multi
	gene2include, patient2include = Multi.white_and_blacklisting(None, None,
		                            None, GENE_BLACKLIST)

	mutation_data = Multi.load_mutation_data_w_cutoff(MUTATION_MATRIX,
		patient2include, gene2include, CUTOFF)
	m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

	# Permute the mutation data
	G  = P.construct_mutation_graph(mutation2patients, patient2mutations)
	Hs = [ P.permute_mutation_data(G.copy(), genes, patients)
	       for i in range(NUM_PERMUTATIONS) ]

	# Calculate the p-value and output the results
	W_prime = sum([ Multi.W(mutation2patients, module, ALPHA)
 		           for module in COLLECTION])
	print  "W' of collection:", W_prime
	print 'P=', matrix_permutation_test(W_prime, Hs, T, K_MIN, K_MAX, ALPHA,
		                                OVERLAPS, SETS_PER_GENE)

if __name__ == "__main__": run()
