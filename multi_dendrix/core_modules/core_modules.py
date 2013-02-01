#!/usr/bin/python

# Load required modules

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Extracts core modules from input collections of gene sets.'\
                  ' Core modules are defined by how often genes appear in the'\
                  ' same gene set together.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--stability_threshold', type=int, default=1,
                        help='Minimum proportion of pathways two genes must '\
                             'both be a member of to be connected in the '\
                             'core modules.')
    parser.add_argument('-i', '--input_collections', nargs="*", required=True,
    	                help='Input files that are collections of gene sets '\
    	                      'output by Multi-Dendrix')
    parser.add_argument('-o', '--output_file', required=True,
    	                help='Output file location.')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args

def load_collection(collection_file):
	"""Extracts the gene sets from a collection file output by Multi-Dendrix

	:type collection_file: string
	:param collection_file: file location of collection to be loaded.

	:returns: **collection** (*list*) - the gene sets (*without* weights) in the input collection file.
	"""
	return [ l.rstrip().split("\t")[1:] for l in open(collection_file) ]

def extract(collections, stability_threshold):
    """Extracts the core modules from a set of collections output by Multi-Dendrix.
    Core modules are defined by how often genes appear in the same gene set together
    for different parameters of Multi-Dendrix ("how often" is tunable using the
    stability threshold).
	
    :type collections: list of lists of lists
    :param collections: multiple collections of gene sets identified by Multi-Dendrix.
    :type stability_threshold: int
    :param stability_threshold: the number of times two genes must appear in the same gene set in order to be grouped in the same module.

    **Returns:**
      A tuple containing the following:
        * **core_modules** (*list of lists*) - the core modules identified by the algorithm
        * **module_graph** (*NetworkX Graph*) - the graph created to identify the core modules, where genes are nodes and pairs of genes are connected by edges weighted by how often they appeared in the same gene set.

    **Example:**
      A simple example with the resulting module graph (edges with weight zero are not shown):
        >>> import networkx as nx
        >>> collections = [ [["G1", "G2", "G3"], ["G4", "G5", "G6"]],
                            [["G1", "G2", "G3"], ["G4", "G5", "G6", "G7"]],
                            [["G1", "G2", "G3"], ["G4", "G5", "G7", "G8"]] ]
        >>> core_modules, module_graph = extract(collections, 1)
        >>> core_modules
        [['G7', 'G6', 'G5', 'G4'], ['G3', 'G2', 'G1']]
        >>> nx.draw_circular(module_graph, node_size=125, font_size=8)

          .. image:: /_static/module_graph.*

    """
    # Load required modules
    import networkx as nx
    
    # Make a set of all genes in any gene set
    genespace = set( [g for P in collections for p in P for g in p ])

    # Determine pathway membership for all pathways
    membership, num_sets = {}, float(len(collections))
    gene2count = dict([(g, 0) for g in genespace])
    for g in genespace:
        membership[g] = dict([(g2, 0) for g2 in genespace if g2 != g])
        gene2count[g] = len([p for P in collections for p in P if g in p])
        
        for g2 in genespace:
            if g2 == g: continue
            for collection in collections:
                for p in collection:
                    if g in p and g2 in p:
                        membership[g][g2] += 1

    # Create a graph using membership proportions as edge weights
    G = nx.Graph()
    G.add_nodes_from(genespace)
    edges = [(g, g2, membership[g][g2])\
             for g in genespace for g2 in membership[g].keys()
             if membership[g][g2] >= stability_threshold]
    G.add_edges_from([(u, v, dict(weight=w)) for u, v, w in edges])
    
    # Return the graph and its connected components
    return [ subG for subG in nx.connected_components(G) if len(subG) > 1 ], G

def run(args):
	"""Identifies the core modules for the given collections, and outputs them to file."""
	# Load collections
	collections = [ load_collection(f) for f in args.input_collections ]
	core_modules, module_graph = extract(collections, args.stability_threshold)

	# Output core modules in single file
	open(args.output_file, 'w').write("\n".join(["\t".join(M)
		                                         for M in core_modules]))

if __name__ == "__main__": run(parse_args())