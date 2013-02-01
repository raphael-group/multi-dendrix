#!/usr/bin/python

# Import globally required modules
import sys
try: import networkx as nx
except ImportError:
	print 'Error!'
	print '\tCould not import NetworkX (http://networkx.github.com).'
	print '\tMake sure NetworkX is in your path.'
	sys.exit(1)
from .. import multi_dendrix as Multi
import random

# Parse args
def parse_args(input_list=None):
	# Parse arguments
	import argparse
	class Args: pass
	args = Args()
	description = 'Creates permuted matrices for a given set of Multi-Dendrix'\
	              ' mutation data parameters using the MeMO permutation '\
	              'method. Requires NetworkX.'
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-m', '--mutation_matrix', required=True,
		                help='File name for mutation data.')
	parser.add_argument('-c', '--cutoff', type=int, default=0, 
		                help='Minimum gene mutation frequency.')
	parser.add_argument('-p', '--patient_whitelist', default=None,
		                help='File of patients to be included.')
	parser.add_argument('-bp', '--patient_blacklist', default=None,
		                help='File of patients to be excluded.')
	parser.add_argument('-g', '--gene_whitelist', default=None,
		                help='File of genes to be included.')
	parser.add_argument('-bg', '--gene_blacklist', default=None,
		                help='File of genes to be excluded.')
	parser.add_argument('-o', '--output_dir', required=True,
		                help='Name of output directory.')
	parser.add_argument('-s', '--start_index', default=1, type=int,
		                help='Start index for name of permuted matrices.')
	parser.add_argument('-n', '--num_matrices', type=int, default=100,
		                help='Number of overlaps allowed per pathway.')
	parser.add_argument('-q', '--Q', type=int, default=100,
		                help='Edge swapping parameter.')
	parser.add_argument('-v', '--verbose', default=False, action='store_true',
		                help='Flag verbose mode.')

	# If called from the command line, parse command line args.
	if input_list: parser.parse_args(input_list, namespace=args)
	else: parser.parse_args(namespace=args)
	
	return args

def log(s):
	sys.stdout.write(s)
	sys.stdout.flush()

def bipartite_double_edge_swap(G, genes, patients, nswap=1, max_tries=1e75):
	"""A modified version of the double_edge_swap function in NetworkX to preserve the bipartite structure of the graph.

	For more details on this function, please see the `original NetworkX function <http://goo.gl/wWxBD>`_ that I shamelessly used to make this one.
	The only major change here is that I ensure that u,v and x,y are of the same type (i.e. genes or patients).

	:type G: NetworkX Graph
	:param G: Bipartite graph G(V, E) representation of mutation data. Vertices are genes and patients, and edges connect genes mutated in particular patients.
	:type genes: list
	:param genes: genes in the mutation data.
	:type patients: list
	:param patients: patients in the mutation data.
	:type nswap: int
	:param nswap: number of edge swaps to perform (default: 1).
	:type max_tries: int
	:param max_tries: maximum number of attempted edge swaps to perform (default: 1e75).

	:returns: Bipartite graph G (modified in place).	

	**See also:** :func:`permute_mutation_data`. 
	"""
	if nswap>max_tries:
		raise nx.NetworkXError("Number of swaps > number of tries allowed.")
	if len(G) < 4:
		raise nx.NetworkXError("Graph has less than four nodes.")
    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
	n=0
	swapcount=0
	keys,degrees=zip(*G.degree().items()) # keys, degree
	cdf=nx.utils.cumulative_distribution(degrees)  # cdf of degree
	while swapcount < nswap:
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
		(ui,xi)=nx.utils.discrete_sequence(2,cdistribution=cdf)
		if ui==xi:
			continue # same source, skip
		u=keys[ui] # convert index to label
		x=keys[xi]
		if (u in genes and x in genes) or (u in patients and x in patients):
			continue # both are genes, skip

		patient1 = u if u in patients else x
		gene1    = x if x in genes else u

		# choose target uniformly from neighbors
		patient2=random.choice( list(G[gene1]) )
		gene2=random.choice( list(G[patient1]) )
	
		# don't create parallel edges
		if (gene1 not in G[patient1]) and (gene2 not in G[patient2]): 
			G.add_edge(gene1,patient1)
			G.add_edge(gene2,patient2)
			
			G.remove_edge(gene1,patient2)
			G.remove_edge(patient1, gene2)
			swapcount+=1
		if n >= max_tries:
			e=('Maximum number of swap attempts (%s) exceeded '%n +
			'before desired swaps achieved (%s).'%nswap)
			raise nx.NetworkXAlgorithmError(e)
		n+=1
	return G

def construct_mutation_graph(mutation2patients, patient2mutations):
	"""Converts mutation data stored as dictionaries into a bipartite NetworkX graph.

    :type mutation2patients: dictionary
    :param mutation2patients: Mapping of genes to the patients in which they are mutated.
    :type patient2mutations: dictionary
    :param patient2mutations: Mapping of patients to the genes they have mutated.

    For more information on the internal format of mutation data used by Multi-Dendrix, see :func:`multi_dendrix.load_mutation_data`.

    :returns: Bipartite NetworkX graph G=(V, E) with genes and patients as nodes, and edges representing a mutation in a particular gene in a particular patient.

    **Examples:**
      A view of example input:
        >>> import networkx as nx
        >>> mutation2patients = {"G1" : ["TCGA-01", "TCGA-02", "TCGA-03"], "G2" : ["TCGA-02"]}
        >>> patient2mutations = {"TCGA-01" : ["G1"], "TCGA-02" : ["G1", "G2"], "TCGA-03" : ["G1"]}
      Simple example of converting mutation data into a bipartite graph:
        >>> G = construct_mutation_graph(mutation2patients, patient2mutations)
        >>> nx.draw_spectral(G)

        .. image:: /_static/mutation_graph.*

    **See also:** :func:`graph_to_mutation_data`, :func:`permute_mutation_data`.

	"""
	genes, patients = mutation2patients.keys(), patient2mutations.keys()
	nodes = genes + patients
	edges = [ (gene, patient) for gene in genes
	          for patient in mutation2patients[gene] ]
	G = nx.Graph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)
	return G

def graph_to_mutation_data(H, genes, patients):
	"""Converts a bipartite NetworkX graph representing mutations in genes in different patients into the mutation data format used by Multi-Dendrix.

	For more information on the mutation data format used by Multi-Dendrix, see :func:`multi_dendrix.load_mutation_data`.

	:type H: NetworkX graph
	:param H: Bipartite graph H(V, E) representation of mutation data. Vertices are genes and patients, and edges connect genes mutated in particular patients.
	:type genes: list
	:param genes: genes in the mutation data.
	:type patients: list
	:param patients: patients in the mutation data.

	:returns: Mutation data tuple in the same format as :func:`multi_dendrix.load_mutation_data`.

	**Examples:**
		A view of example input:
			>>> import networkx as nx
			>>> H = nx.Graph()
			>>> H.add_edges_from([("G1", "TCGA-01"), ("G1", "TCGA-02"), ("G1", "TCGA-03"),
	    		                  ("G2", "TCGA-02")])
			>>> nx.draw_spectral(H)

			.. image:: /_static/mutation_graph.*
		Converting the graph into Multi-Dendrix mutation data format:
			>>> graph_to_mutation_data(H, ["G1", "G2"], ["TCGA-01", "TCGA-02", "TCGA-03"])
			>>> (2, 3, ['G2', 'G1'], ['TCGA-03', 'TCGA-02', 'TCGA-01'],
				{'G2': set(['TCGA-02']), 'G1': set(['TCGA-03', 'TCGA-02', 'TCGA-01'])},
				{'TCGA-03': set(['G1']), 'TCGA-02': set(['G2', 'G1']), 'TCGA-01': set(['G1'])})

	**See also:** :func:`permute_mutation_data`, :func:`construct_mutation_graph`.

	"""
	mutation2patients, patient2mutations = dict([(g, set()) for g in genes]), dict( )
	for patient in patients:
		mutations = H[patient]
		patient2mutations[patient] = set( mutations )
		for g in mutations: mutation2patients[g].add( patient )

	genes, patients = mutation2patients.keys(), patient2mutations.keys()
	m, n = len(genes), len(patients)
	return m, n, genes, patients, mutation2patients, patient2mutations

def permute_mutation_data(G, genes, patients, Q=100):
	"""Permutes the given mutation data stored in bipartite graph G=(V, E) by performing | E | * Q edge swaps.

	:type G: NetworkX Graph
	:param G: Bipartite graph G=(V, E) representation of mutation data. Vertices are genes and patients, and edges connect genes mutated in particular patients.
	:type genes: list
	:param genes: genes in the mutation data.
	:type patients: list
	:param patients: patients in the mutation data.
	:type Q: int
	:param Q: constant multiplier for number Q * | E | of edge swaps to perform (default and suggested value: 100). See `Milo et al. (2003) <http://goo.gl/d723i>`_ for details on choosing Q.

	:returns: Permuted version of G reformatted  into the mutation data format used by Multi-Dendrix (see :func:`graph_to_mutation_data` and :func:`multi_dendrix.load_mutation_data`).

	**Examples:**
	  A view of example input:
	    >>> import networkx as nx
	    >>> G = nx.Graph()
	    >>> G.add_edges_from([("G1", "TCGA-01"), ("G1", "TCGA-02"), ("G1", "TCGA-03"),
	    	("G2", "TCGA-02"), ("G3", "TCGA-01"), ("G3", "TCGA-02"), ("G4", "TCGA-03")])
	    >>> nx.draw_spectral(G, dpi=72, node_size=125, font_size=8)
	    
	    .. image:: /_static/permute_mutation_data_before.*
	  Permute the mutation data:
	    >>> M = permute_mutation_data(G, ["G1", "G2", "G3", "G4"], ["TCGA-01", "TCGA-02", "TCGA-03"])
	    >>> M
	    (4, 3, ['G4', 'G3', 'G2', 'G1'], ['TCGA-03', 'TCGA-02', 'TCGA-01'],
	    	{'G4': set(['TCGA-02']), 'G3': set(['TCGA-02', 'TCGA-01']), 'G2': set(['TCGA-03']),
	    	'G1': set(['TCGA-03', 'TCGA-02', 'TCGA-01'])},
	    	{'TCGA-03': set(['G2', 'G1']), 'TCGA-02': set(['G4', 'G3', 'G1']),
	    	'TCGA-01': set(['G3', 'G1'])})
	    >>> H = construct_mutation_graph(M[-1], M[-2])
	    >>> nx.draw_spectral(H, dpi=72, node_size=125, font_size=8)

	    .. image:: /_static/permute_mutation_data_after.*
	**See also:** :func:`construct_mutation_graph`, :func:`graph_to_mutation_data`.
	"""

	H = G.copy()
	bipartite_double_edge_swap(H, genes, patients, nswap=Q * len( G.edges() ))
	return graph_to_mutation_data(H, genes, patients)

def run(args):
	"""Permutes the given mutation data a given number of times."""
	# Load mutation data using Multi-Dendrix and output as a temporary file
	if args.verbose: log('Loading mutation data...')
	
	include = Multi.white_and_blacklisting(args.patient_whitelist,
		args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
	gene2include, patient2include = include

	mutation_data = Multi.load_mutation_data_w_cutoff(args.mutation_matrix,
		patient2include, gene2include, args.cutoff)
	m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

	if args.verbose: log('done!\n\n')

	# Make sure output directory exists
	import os
	os.system('mkdir -p ' + args.output_dir)

	# Construct bipartite graph from mutation data
	if args.verbose: log('Creating bipartite graph...')
	
	G = construct_mutation_graph(mutation2patients, patient2mutations)

	if args.verbose:
		log('done!\n\n')
		print 'Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.\n'

	# Create permuted matrices and save to file
	for i in range(args.num_matrices):
		if args.verbose: log('+')

		# Permute bipartite graph and output as a patient adjacency list
		mutation_data = permute_mutation_data(G, genes, patients, args.Q)
		_, _, _, _, mutation2patients, patient2mutations = mutation_data
		adj_list = [ p + "\t" + "\t".join( patient2mutations[p] )
		             for p in patients ]

		filename = args.output_dir + "/" + str(i + args.start_index) + '.txt'
		open(filename, 'w').write('\n'.join(adj_list))

	if args.verbose: print

if __name__ == "__main__": run(parse_args())
