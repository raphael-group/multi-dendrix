#!/usr/bin/python
"""
Direct interactions test
============================
"""

def run():
	"""
	This script performs the direct interaction test on a collection of gene 
	sets. It gives an example of how one could use various functions from 
	different modules in order to evaluate the collection without obtaining
	the collection from the Multi-Dendrix pipeline. The program first loads
	the iRef network, permutes it the given number of times, and then performs
	the direct interactions test of the pathway set.

	The default parameters permute the iRefIndex *once*, and only swaps | E |
	total edges (i.e. Q=1). On my machine, this script runs in 19.2 seconds.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../'))
	import multi_dendrix.permute.ppi_network as P
	import multi_dendrix.evaluate.network as Test

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	COLLECTION = [["CDKN2A", "RB1", "CDK4"], ["TP53", "MDM2", "MDM4"]]
	NETWORK_EDGELIST = "../network/iref_edgelist"
	NUM_PERMUTATIONS=1
	Q=1
	USE_DISTANCE=False # Use average pairwise distance as the test statistic

	# Load the network and permute it
	G  = P.load_network(NETWORK_EDGELIST)
	Hs = [ P.permute_network(G, Q) for i in range(NUM_PERMUTATIONS) ]

	# Perform direct interactions test
	results = Test.evaluate_collection(COLLECTION, G, Hs, USE_DISTANCE)
	test_name, set_statistic, set_pval, module_results = results

	# Output pathway set evaluation
	print 'Collection test'
	print '\t' + test_name + ':', set_statistic
	print '\tP-value:', set_pval
	
	# Output individual pathway evaluation
	print '\nIndividual gene sets test'
	for p, stat, pval in module_results:
		print '\tGenes:', ', '.join(p)
		print '\tStatistic:', stat
		print '\tP-value:', pval, '\n'



if __name__ == "__main__": run()
