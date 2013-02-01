#!/usr/bin/python

# Load required modules
import sys, os, networkx as nx
from itertools import combinations, product

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Calculates network permutation test for given collection '\
                  'of gene sets on the given network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-n', '--network_edgelist', required=True,
                        help='PPI edgelist location.')
    parser.add_argument('-i', '--permuted_networks_dir', required=True,
                        help='Directory of permuted networks.')
    parser.add_argument('-p', '--collection_file', required=True,
                        help='File containing input collection of gene sets.')
    parser.add_argument('-o', '--output_file', default=None,
                        help='Name of output file.')
    parser.add_argument('-d', '--distance', default=False, action='store_true',
                        help='Flag average pairwise distance test.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Flag verbose mode.')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args

###############################################################################
# Helper functions for performing the direct interactions test

def load_collection(collection_file):
	"""Extracts the gene sets from a collection file output by Multi-Dendrix.

	:type collection_file: string
	:param collection_file: file location of collection to be loaded.

	:returns: **collection** (*list*) - the gene sets (*without* weights) in the input collection file.
	"""
	return [ l.rstrip().split("\t")[1:] for l in open(collection_file) ]

# Wrapper for NetworkX's read_edgelist function
def load_network(network_edgelist): return nx.read_edgelist(network_edgelist)

def remove_name_annotation(genes):
	"""Removes annotation of genes or mutation classes for CNAs.
	By convention, I annotated CNA mutation classes with '(A)' for amplification
	and '(D)' for deletion. This removes these annoations so the genes can be mapped
	to PPI networks. For more information on naming conventions of genes and 
	mutation classes, see :doc:`/file_formats`.

	:type genes: list
	:param genes: gene names or IDs.

	:returns: **genes** (*list*) - list of genes with annotation removed.
	"""
	return map(lambda g: g.replace("(A)", "").replace("(D)", ""), genes)

# Wrapper for itertools combos function
def combos(xs, n): return list(combinations(xs, n))

# Wrapper for itertools product function
def prod(xs, ys): return list(product(xs, ys))

def pairs_within(collection):
    """Given a collection of gene sets, returns all :math:`{n \choose 2}` pairs from within EACH gene set."""
    return reduce(lambda pairs, p: pairs + combos(p, 2), collection, [])

def pairs_between(collection):
    """Given a collection of gene sets, returns all pairs from each distinct pair of gene sets."""
    pairs_of_gsets = combos(collection, 2)
    return reduce(lambda pairs, (p1, p2): pairs + prod(p1, p2), pairs_of_gsets, [])

###############################################################################
# Functions for performing the average pairwise distance test
def avg_pair_dist_of_gene_set(network, gset):
	"""Counts the number of interactions among genes in the given set.

	:type network: NetworkX Graph
	:param network: input PPI network.
	:type gset: list
	:param gset: genes to test.

	:returns: **avg_pair_dist** (*int*) - average pairwise distance among genes in *gset*.

	**See also:** :func:`eval_gene_sets_by_dist`, :func:`sum_dist`.

	**Note:** Used by :func:`eval_gene_sets_by_dist` when performing the *average pairwise distance test* on individual gene sets.
	"""	
	min_dist = 1e100
	gset = list( set( gset ) )
	k = len( gset )
	pairs = [ (gset[i], gset[j]) for i in range(k) for j in range(i+1, k) ]
	return float( sum_dist(network, pairs) ) / float( len(pairs ) )

def eval_gene_sets_by_dist(collection, G, Hs):
	"""Evaluate each of the *t* individual gene sets by comparing the number of interactions within a gene set in the original PPI network G to the permuted networks Hs.

	:type collection: list of lists
	:param collection: Collection of gene sets to be tested.
	:type G: NetworkX Graph.
	:param G: Original PPI network.
	:type Hs: list of NetworkX Graphs
	:param Hs: permuted versions of G.

	**Returns:**
	  A list of *t* tuples where each tuple contains:
	    * **gset** (*list*) - the original gene set.
	    * **avg_pair_dist** (*float*) - the average pairwise  genes set in the original network.
	    * **pval**  (*float*) - the empirical *p*-value of the gene set's test statistic (i.e. the average pairwise distance) in the permuted networks.

	**See also:** :func:`evaluate_collection`, :func:`eval_gene_sets_by_interactions`, :func:`avg_pair_dist_test`.
	"""
	stats, pvals = [], []
	for gset in collection:
		avg_pair_dist = avg_pair_dist_of_gene_set(G, gset)
		permuted_dists  = [ avg_pair_dist_of_gene_set(H, gset) for H in Hs ]
		extreme = [ n for n in permuted_dists if n <= avg_pair_dist]
		pval = float( len(extreme) ) / float( len(Hs) )

		stats.append( avg_pair_dist )
		pvals.append( pval )

	return zip(collection, stats, pvals)

def sum_dist(network, pairs):
	"""Given a PPI and a list of gene pairs, returns the sum of the shortest paths between each pair.
	
	:type network: NetworkX Graph
	:param network: input PPI network.
	:type pairs: list of tuples
	:param pairs: pairs of genes to test.

	:returns: **distance** (*int*) - the sum of the distances between each gene pair.

	**See also:** :func:`avg_pair_dist_ratio`, :func:`dist`.
	"""
	return sum( [ dist(network, g1, g2) for g1, g2 in pairs ] )

def dist(network, g1, g2):
	"""Returns the length of the shortest path between g1 and g2 in PPI. If no path exists, returns 1e100.

	:type network: NetworkX Graph
	:param network: input PPI network.
	:type pairs: list of tuples
	:param pairs: pairs of genes to test.

	:returns: **distance** (*float*) - distance between the two genes on the network. If no path between the genes exists, returns 1e100.

	**See also:** :func:`direct_interactions_stat`, :func:`interact`.

    """
	try: return len(nx.shortest_path(network, g1, g2))-1
	except (nx.NetworkXNoPath, nx.NetworkXError): return 1e100

def avg_pair_dist_ratio(network, collection):
	"""Calculates the ratio of the average pairwise distance among genes within the *same* gene set and genes within *different* gene sets.

	:type network: NetworkX Graph.
	:param network: PPI network (original or permuted).
	:type: collection: list of lists
	:param collection: collection of gene sets to be tested.

	:returns: :math:`\\rho(\mathbf{P})` (*float*) - see definition below.

	**Test statistic**:
	  For a collection of gene sets :math:`\mathbf{P}`:
	    * Let :math:`\\alpha'(\mathbf{P})=` average pairwise distance among pairs of genes in the *same* gene set.
	    * Let :math:`\\beta'(\mathbf{P})=` average pairwise distance among pairs of genes in *different* gene sets.
	  Then the average pairwise distance ratio :math:`\\rho(\mathbf{P}) = \\alpha'(\mathbf{P})-\\beta'(\mathbf{P})`. *Lower* values of :math:`\\rho(\mathbf{P})` indicate stronger collections.

	  *Note*: if no path exists between any genes in the collection, the average pairwise distance ratio will go to infinity since the distance between two nodes in different components (or nodes not in the graph) is ill-defined.

	**Examples:**
	  A view of example input:
	    >>> import networkx as nx
	    >>> collection = [["G1", "G2"], ["G3", "G4"]]
	    >>> G = nx.Graph()
	    >>> G.add_edges_from([("G1", "G2"), ("G3", "G4"), ("G2", "G5"), ("G4", "G5")])
	    >>> nx.draw_circular(G, node_size=125, font_size=8)

	    .. image:: /_static/direct_interactions_test_graph.png
	  Calculating :math:`\\alpha'(\mathbf{P})` and :math:`\\beta'(\mathbf{P})`:
	    >>> dist_within = sum_dist(G, pairs_within(collection))
	    >>> alpha_prime  = float(dist_within) / float(len(pairs_within(collection)))
	    >>> (dist_within, alpha_prime)
	    (2, 1.0)
	    >>> dist_between = sum_dist(G, pairs_between(collection))
	    >>> beta_prime  = float(dist_between) / float(len(pairs_between(collection)))
	    >>> (dist_between, beta_prime)
	    (12, 3.0)
	  The direct interactions test statistic:
	    >>> avg_pair_dist_ratio(G, collection)
	    0.33

	**See also:** :func:`direct_interactions_test`, :func:`eval_gene_sets_by_interactions`.
	"""
	# Ensure all genes are in the network
	if not all([ g in network.nodes() for P in collection for g in P ]):
		return 1e100

	# Enumerate pairs of genes within the same and between gene sets
	within, between = pairs_within(collection), pairs_between(collection)
	within_dists    = [ dist(network, g1, g2) for g1, g2 in within ]
	between_dists   = [ dist(network, g1, g2) for g1, g2 in between ]

	# Ensure all genes in the same gene set are in the same component
	if not all([ d != 1e100 for d in within_dists ]):
		return 1e100

	# Check if any genes in different gene sets are not in the same component
	if not all([ d != 1e100 for d in between_dists ]):
		return 0

	# If all checks out, return the ratio
	d_within  = float(sum(within_dists)) / float(len(within))
	d_between = float(sum(between_dists)) / float(len(between))
	return d_within / d_between

def avg_pair_dist_test(collection, G, Hs):
	"""Performs the *average pairwise distance test* on the collection of gene sets.

	:type collection: list of lists
	:param collection: collection of gene sets to be tested.
	:type G: NetworkX Graph.
	:param G: original PPI network.
	:type Hs: list of NetworkX Graphs
	:param Hs: permuted versions of G.

	**Returns:**
	  A tuple consisting of the following:
	    * :math:`\\rho(\mathbf{P})` (*float*) - the average pairwise distance test statistic.
	    * **pval** (*float*) - the average pairwise distance test *p*-value of the given collection.

	**Test:**

	For the definition of the test statistic used by the *average pairwise distance test*, see :func:`avg_pair_dist_ratio`.
	The *average pairwise distance test* compares the average pairwise distance test statistic of the collection on the original and permuted PPI networks.

	**Note:** Genes that are not in the graph *G* or are in different connected components than the other genes are given an average pairwise distance of 1e100. These genes should be **removed** before using the *average pairwise distance test*.

	**Examples:**
	  A view of example input:
	    >>> import networkx as nx
	    >>> collection = [["G1", "G2"], ["G3", "G4"]]
	    >>> G = nx.Graph()
	    >>> G.add_edges_from([("G1", "G2"), ("G3", "G4"), ("G2", "G5"), ("G4", "G5")])
	    >>> Hs =[ nx.double_edge_swap(G.copy(), nswap=10, max_tries=1e75) for i in range(10) ]
	    >>> nx.draw_circular(G, node_size=125, font_size=8)
	  
	    .. image:: /_static/direct_interactions_test_graph.png
	  A simple example
	    >>> avg_pair_dist_test(collection, G, Hs)
	    (0.33, 0.10)

	**See also:** :func:`avg_pair_dist_ratio`, :func:`eval_gene_sets_by_dist`, :func:`direct_interactions_test`.
	"""
	# Remove name annotation for genes in the gene sets
	collection = [ list(set(remove_name_annotation(P))) for P in collection ]
	
	# Calculate the average pairwise distance ratio in the original network
	# and the permuted networks. Lower ratios indicate stronger collections.
	ratio  = avg_pair_dist_ratio(G, collection)
	permuted_ratios = [ avg_pair_dist_ratio(H, collection) for H in Hs ]
	count = len( [r for r in permuted_ratios if r <= ratio ] )
	pval  = float( count ) / float( len( Hs ) )

	return ratio, pval

###############################################################################
# Functions for performing the direct interactions test

def num_interactions_in_gene_set(network, gset):
	"""Counts the number of interactions among genes in the given set.

	:type network: NetworkX Graph
	:param network: input PPI network.
	:type gset: list
	:param gset: genes to test.

	:returns: **count** (*int*) - number of interactions among genes in *gset*.

	**See also:** :func:`eval_gene_sets_by_interactions`, :func:`count_interactions`.

	**Note:** Used by :func:`eval_gene_sets_by_interactions` when performing the *direct interactions test* on individual gene sets.
	"""
	best_count = 0
	gset = list( set( gset ) )
	k = len( gset )
	pairs = [ (gset[i], gset[j]) for i in range(k) for j in range(i+1, k) ]
	return count_interactions(network, pairs)

def eval_gene_sets_by_interactions(collection, G, Hs):
	"""Evaluate each of the *t* individual gene sets by comparing the number of interactions within a gene set in the original PPI network G to the permuted networks Hs.

	:type collection: list of lists
	:param collection: collection of gene sets to be tested.
	:type G: NetworkX Graph.
	:param G: original PPI network.
	:type Hs: list of NetworkX Graphs
	:param Hs: permuted versions of G.

	**Returns:**
	  A list of *t* tuples where each tuple contains:
	    * **gset** (*list*) - the original gene set.
	    * **num_interactions** (*int*) - the number of interactions among genes set in the original network.
	    * **pval**  (*float*) - the empirical *p*-value of the gene set's test statistic (i.e. the number of interactions) in the permuted networks.

	**See also:** :func:`evaluate_collection`, :func:`eval_gene_sets_by_dist`, :func:`direct_interactions_test`.
	"""
	stats, pvals = [], []
	for gset in collection:
		num_interactions = num_interactions_in_gene_set(G, gset)
		permuted_counts  = [ num_interactions_in_gene_set(H, gset) for H in Hs]
		extreme = [ n for n in permuted_counts if n >= num_interactions]
		pval = float( len(extreme) ) / float( len(Hs) )

		stats.append( num_interactions )
		pvals.append( pval )

	return zip(collection, stats, pvals)

def interact(network, g1, g2):
	"""Returns true if g1 interacts with g2 in PPI.

	:type network: NetworkX Graph
	:param network: input PPI network.
	:type g1: string
	:param g1: gene in PPI.
	:type g2: string
	:param g2: gene in PPI.

	:returns: True if g1 interactions with g2 in the given network, False otherwise.
	"""
	try: return g2 in nx.neighbors(network, g1)
	except nx.exception.NetworkXError: return False

def count_interactions(network, pairs):
	"""Given a PPI and a list of gene pairs, returns the number of genes that interact.
	
	:type network: NetworkX Graph
	:param network: input PPI network.
	:type pairs: list of tuples
	:param pairs: pairs of genes to test.

	:returns: **count** (*int*) - the number of gene pairs that interact.

	**See also:** :func:`direct_interactions_stat`, :func:`interact`.
	"""
	return sum( [int(interact(network, g1, g2)) for g1, g2 in pairs])

def direct_interactions_stat(network, collection):
	"""Calculates the difference of the normalized number of interactions among genes within the *same* gene set and genes within *different* gene sets.

	:type network: NetworkX Graph.
	:param network: PPI network (original or permuted).
	:type: collection: list of lists
	:param collection: collection of gene sets to be tested.

	:returns: :math:`\\nu(\mathbf{P})` (*float*) - see definition below.

	**Test statistic**:
	  For a collection of gene sets :math:`\mathbf{P}`:
	    * Let :math:`\\alpha(\mathbf{P})=` normalized number of interactions among pairs of genes in the *same* gene set.
	    * Let :math:`\\beta(\mathbf{P})=` normalized number of interactions among pairs of genes in *different* gene sets.
	  Then the direct interactions statistic :math:`\\nu(\mathbf{P}) = \\alpha(\mathbf{P})-\\beta(\mathbf{P})`. Higher values of :math:`\\nu(\mathbf{P})` indicate stronger collections.

	**Examples:**
	  A view of example input:
	    >>> import networkx as nx
	    >>> collection = [["G1", "G2"], ["G3", "G4"]]
	    >>> G = nx.Graph()
	    >>> G.add_edges_from([("G1", "G2"), ("G3", "G4"), ("G2", "G5"), ("G4", "G5")])
	    >>> nx.draw_circular(G, node_size=125, font_size=8)

	    .. image:: /_static/direct_interactions_test_graph.png
	  Calculating :math:`\\alpha(\mathbf{P})` and :math:`\\beta(\mathbf{P})`:
	    >>> int_within = count_interactions(G, pairs_within(collection))
	    >>> alpha  = float(int_within) / float(len(pairs_within(collection)))
	    >>> (int_within, alpha)
	    (2, 1.0)
	    >>> int_between = count_interactions(G, pairs_between(collection))
	    >>> beta  = float(int_between) / float(len(pairs_between(collection)))
	    >>> (int_between, beta)
	    (0, 0)
	  The direct interactions test statistic:
	    >>> direct_interactions_stat(G, collection)
	    1.0

	**See also:** :func:`direct_interactions_test`, :func:`eval_gene_sets_by_interactions`.
	"""
	within_pairs = pairs_within(collection)
	interactions_within  = count_interactions(network, within_pairs)
	alpha = float(interactions_within) / float(len(within_pairs))

	between_pairs = pairs_between(collection)
	interactions_between = count_interactions(network, between_pairs)
	beta = float( interactions_between ) / float(len(between_pairs))

	return alpha - beta

def direct_interactions_test(collection, G, Hs):
	"""Performs the *direct interactions test* on the collection of gene sets.

	:type collection: list of lists
	:param collection: collection of gene sets to be tested.
	:type G: NetworkX Graph.
	:param G: original PPI network.
	:type Hs: list of NetworkX Graphs
	:param Hs: permuted versions of G.

	**Returns:**
	  A tuple consisting of the following:
	    * :math:`\\nu(\mathbf{P})` (*float*) - the direct interactions test statistic.
	    * **pval** (*float*) - the direct interactions test *p*-value of the given collection.

	**Test:**

	For the definition of the test statistic used by the *direct interactions test*, see :func:`direct_interactions_stat`.
	The *direct interactions test* compares the direct interactions test statistic of the collection on the original and permuted PPI networks.

	**Examples:**
	  A view of example input:
	    >>> import networkx as nx
	    >>> collection = [["G1", "G2"], ["G3", "G4"]]
	    >>> G = nx.Graph()
	    >>> G.add_edges_from([("G1", "G2"), ("G3", "G4"), ("G2", "G5"), ("G4", "G5")])
	    >>> Hs =[ nx.double_edge_swap(G.copy(), nswap=10, max_tries=1e75) for i in range(10) ]
	    >>> nx.draw_circular(G, node_size=125, font_size=8)
	  
	    .. image:: /_static/direct_interactions_test_graph.png
	  A simple example
	    >>> direct_interactions_test(collection, G, Hs)
	    (1.0, 0.1)

	**See also:** :func:`direct_interactions_stat`, :func:`eval_gene_sets_by_interactions`, :func:`avg_pair_dist_test`.
	"""
	# Remove name annotation for genes in the gene sets
	collection = [ list(set(remove_name_annotation(P))) for P in collection ]
	
	# Calculate the number of interactions statistic in the original network
	# and the permuted networks. Higher statistics indicate stronger
	# collections.
	stat  = direct_interactions_stat(G, collection)
	permuted_stats = [ direct_interactions_stat(H, collection) for H in Hs ]
	count = len( [s for s in permuted_stats if s >= stat ] )
	pval  = float( count ) / float( len( Hs ) )

	return stat, pval

###############################################################################
# Main functions for evaluating collections of gene sets

def evaluate_collection(collection, G, Hs, distance=False):
	"""Given a collection of gene sets, the original network, and a set of permuted networks, calculates the empirical *p*-value using the direct interactions statistic (or average pairwise distance statistic).
	For details on the tests performed on collections of gene sets (with examples), see :func:`direct_interactions_test` (or :func:`avg_pair_dist_test`).

	:type collection: list of lists
	:param collection: collection of gene sets to be tested.
	:type G: NetworkX Graph.
	:param G: original PPI network.
	:type Hs: list of NetworkX Graphs
	:param Hs: permuted versions of G.
	:type distance: bool
	:param distance: Flag whether to use average pairwise distance as the test statistic.

	**Returns:**
	  * **test_name** (*string*) - name of test (either "Direct_Interactions_Statistic" or "Average_Pairwise_Distance_Ratio").
	  * **statistic** (*float*) - test statistic of the collection on the original network
	  * **pval** (float) - empirical *p*-value of the collection on the original versus permuted networks.
	  * **gene_set_results** (list) - evaluation of each gene set using the original and permuted networks (see :func:`eval_gene_sets_by_interactions` and :func:`eval_gene_sets_by_dist` for details).

	**See also:** :func:`direct_interactions_test`, :func:`avg_pair_dist_test`.
	"""
	# Calculate test statistic and p-value
	if distance:
		test_name = 'Average_Pairwise_Distance_Ratio'
		statistic, pval = avg_pair_dist_test(collection, G, Hs)
		gene_set_results = eval_gene_sets_by_dist(collection, G, Hs)
	else:
		test_name = 'Direct_Interactions_Statistic'
		statistic, pval = direct_interactions_test(collection, G, Hs)
		gene_set_results = eval_gene_sets_by_interactions(collection, G, Hs)

	return test_name, statistic, pval, gene_set_results

def run(args):
	"""Performs the direct interactions test (or average pairwise distance test) on a the input collection of gene sets using the empirical distribution of input permuted PPI networks."""
	# Load input
	if args.verbose: print 'Loading networks...'
	collection = load_collection(args.collection_file)
	G         = load_network(args.network_edgelist)
	Hs        = [ load_network(args.permuted_networks_dir + "/" + fh)
	              for fh in os.listdir(args.permuted_networks_dir) ]
	
	# Evaluate the collection
	if args.verbose: print 'Evaluating input collection...'
	results = evaluate_collection(collection, G, Hs, args.distance)
	test_name, statistic, pval, gene_set_results = results	

	# Output and return results
	if args.output_file:
		header = 'Testing\t' + test_name + '\tP-Value\n'
		output = [ '\t'.join(['Collection', str(statistic), str(pval)]) ]
		for gene_set, stat, p in gene_set_results:
			output.append( '\t'.join([', '.join(gene_set), str(stat), str(p)]) )
		open(args.output_file, 'w').write(header + '\n'.join(output))
	
	if args.verbose:
		print 'Collection of gene sets test'
		print '\t' + test_name + ':', statistic
		print '\tP-value:', pval
		print '\nIndividual gene sets test'
		for p, stat, pval in gene_set_results:
			print '\tGenes:', ', '.join(p)
			print '\tStatistic:', stat
			print '\tP-value:', pval, '\n'

	return statistic, pval, gene_set_results

if __name__ == "__main__": run(parse_args())