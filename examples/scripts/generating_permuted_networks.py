#!/usr/bin/python
"""
Generating permuted networks
============================
"""

def run():
	"""
	This is a simple Python script for generating permuted PPI networks. The
	default parameter values are for permuting iRefIndex PPI network in
	examples/networks/ directory, which is a large PPI network.
	The main utility of this script is that it is much faster to generate
	permutation data before running the permutation test.

	On my machine, it takes 14 minutes and 5 seconds to permute iRef **once**
	using Q=100, since this requires performing 21,274,600 edge swaps.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../'))
	import multi_dendrix.permute.ppi_network as P

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	NETWORK_EDGELIST="../network/iref_edgelist"
	OUTPUT_DIR="./sample_output/permuted_networks"
	NUM_MATRICES=1

	args = [ "-n", NETWORK_EDGELIST, "-o", OUTPUT_DIR, "-v" ]
	if NUM_MATRICES: args += [ "-c", NUM_MATRICES ]

	# Use the permute_mutation_data module to generate the permuted matrices
	P.run(P.parse_args(map(str, args)))

if __name__ == "__main__": run()
