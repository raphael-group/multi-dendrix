#!/usr/bin/python

# Load required modules
import networkx as nx

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Creates permuted networks for the given network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-n', '--network_edgelist', required=True,
                        help='PPI edgelist location.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Name of output directory.')
    parser.add_argument('-q', '--Q', default=100, type=int,
                        help='Edge swap constant.')
    parser.add_argument('-c', '--num_networks', default=100, type=int,
                        help='Number of permuted networks to create.')
    parser.add_argument('-s', '--start_index', default=1, type=int,
                        help='Start index for name of permuted networks.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Flag verbose mode.')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args

def load_network(network_file): 
    """Wrapper for the NetworkX `read_edgelist <http://goo.gl/xfDPS>`_ function."""
    return nx.read_edgelist(network_file)

def permute_network(G, Q):
    """Permutes the given graph G=(V, E) by performing | E | * Q edge swaps.

    :type G: NetworkX Graph
    :param G: PPI network.
    :type Q: int
    :type Q: int
    :param Q: constant multiplier for number Q * | E | of edge swaps to perform (default and suggested value: 100). See `Milo et al. (2003) <http://goo.gl/d723i>`_ for details on choosing Q.

    :returns: a permuted version of H (G is not modified)

    **Examples:**
      A view of example input:
        >>> import networkx as nx
        >>> G = nx.Graph()
        >>> G.add_edges_from([["G1", "G2"], ["G2", "G3"], ["G2", "G4"], ["G1", "G4"],
                              ["G2", "G5"], ["G5", "G6"], ["G5", "G7"]])
        >>> nx.draw_spectral(G, node_size=125, font_size=8)

        .. image:: /_static/ppi_network.*
      Permuting the network by performing | E | * 10 edge swaps.
        >>> permute_network(G, 10)

        .. image:: /_static/permuted_ppi_network.*

    **See also:** :func:`load_network`, :func:`multi_dendrix.permute.mutation_data.permute_mutation_data`.

    **Notes:** Uses the NetworkX `double_edge_swap <http://goo.gl/wWxBD>`_ function.

    """
    H = G.copy()
    nx.double_edge_swap(H, nswap=Q*len( G.edges() ), max_tries=1e75)
    return H

def run(args):
    """Permutes the given PPI network the specified number of times."""
    import sys, os
    # Load network
    G = load_network(args.network_edgelist)

    if args.verbose:
        print 'Input network has', len( G.edges() ), 'edges among', len(G.nodes()),
        print 'nodes.\nPerforming', len( G.edges() ) * args.Q, 'edge swaps.'

    # Make sure output directory exists
    os.system('mkdir -p ' + args.output_dir)

    # Permute network and output files
    for i in range(args.num_networks):
        if args.verbose:
            sys.stdout.write('+')
            sys.stdout.flush()

        # Permute graph and output as an edge list
        H = permute_network(G, args.Q)
        filename = args.output_dir + "/" + str(i + args.start_index) + ".txt"
        nx.write_edgelist(H, filename)

    if args.verbose: print

if __name__ == "__main__": run(parse_args())
