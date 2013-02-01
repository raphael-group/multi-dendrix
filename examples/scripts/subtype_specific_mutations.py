#!/usr/bin/python
"""
Analyzing subtype-specific mutations
======================================
"""

def run():
	"""
	This is a simple Python script for computing *p*-values of subtype-specific
	mutations. The script requires mutation data and a list of patients and
	their respective subtypes, The default parameters here use the BRCA 
	mutation data in the examples/mutation_data/BRCA directory.

	On my machine, it takes 0.5 seconds to compute the subtype-specific
	mutations on the BRCA dataset.
	"""
	# Load required modules
	import sys, os
	sys.path.insert(1, os.path.abspath('../../'))
	import multi_dendrix.subtypes as Sub

	# Default parameters for this tutorial (CHANGE AT YOUR OWN RISK)
	MUTATION_MATRIX="../mutation_data/BRCA/BRCA.m2"
	PATIENT_LIST="../mutation_data/BRCA/BRCA_subtypes.lst"
	GENE_LIST="../mutation_data/BRCA/BRCA.glst"
	OUTPUT_FILE=None

	args = [ "-m", MUTATION_MATRIX, "-p", PATIENT_LIST, "-v" ]
	if GENE_LIST: args += [ "-g", GENE_LIST ]
	if OUTPUT_FILE: args += [ "-o", OUTPUT_FILE ]

	# Use the permute_mutation_data module to generate the permuted matrices
	Sub.run(Sub.parse_args(map(str, args)))

if __name__ == "__main__": run()
