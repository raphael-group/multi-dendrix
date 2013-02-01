#!/usr/bin/python

# Load required modules
from .. import multi_dendrix as Multi

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Calculates network permutation test for given collection '\
                  'of gene sets on the given network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--permuted_matrices_dir', required=True,
                        help='Directory of permuted networks.')
    parser.add_argument('-p', '--collection_file', required=True,
                        help='File containing input collection of gene sets.')
    parser.add_argument('-o', '--output_file', default=None,
                        help='Name of output file.')
    parser.add_argument('-k_min', '--min_gene_set_size', required=True, type=int,
                        help='Minimum gene set size.')
    parser.add_argument('-k_max', '--max_gene_set_size', required=True, type=int,
                        help='Maximum gene set size.')
    parser.add_argument('-t', '--num_pathways', required=True, type=int,
                        help='Number of gene sets.')
    parser.add_argument('-a', '--alpha', type=float, default=1.0,
                        help='Parameter that changes weight function W by '\
                        'weighting the penalty of coverage overlap.')    
    parser.add_argument('--delta', type=int, default=0,
                        help='Number of overlaps allowed per gene set.')    
    parser.add_argument('--lmbda', type=int, default=1,
                        help='Number of gene sets a gene can be a member of.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Flag verbose mode.')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args

def load_permuted_matrices(input_dir):
    """Loads all files as mutation data in a directory of permuted mutation data.

    :type input_dir: string
    :param input_dir: location of mutation data directory.

    :returns: a list of mutation data tuples (for more information on the mutation data format used internally by Multi-Dendrix, see :func:`multi_dendrix.load_mutation_data`).

    **See also:** :func:`matrix_permutation_test`.
    """
    from os import listdir
    return [ Multi.load_mutation_data(input_dir + "/" + fh)
	         for fh in listdir(input_dir) ]

def load_w_prime(collection_file):
    """Loads the weights and sums them from a Multi-Dendrix output file.

    :type collection_file: string
    :param collection_file: file location of collection output by Multi-Dendrix.

    :returns: *W'(M)*, where *M* is the input collection of gene sets.

    **See also:** :func:`matrix_permutation_test`, :func:`load_permuted_matrices`.

    """
    return sum([int(l.rstrip().split("\t")[0]) for l in open(collection_file)])

def matrix_permutation_test(W_prime, permuted_matrices, t, k_min, k_max, alpha,
                           delta, lmbda):
    """Computes the statistical significance of a collection found by Multi-Dendrix using an empricial distribution of mutation data.

    :type W_prime: int
    :param W_prime: sum of the weights of an optimal collection of gene sets identified by Multi-Dendrix. This function compares this score to the optimal score on permuted mutation data.
    :type permuted_matrices: list
    :param permuted_matrices: contains the mutation data tuples loaded by :func:`load_permuted_matrices`.
    :type t: int
    :param t: number of gene sets
    :type k_min: int
    :param k_min: minimum gene set size.
    :type k_max: int
    :param k_max: maximum gene set size.
    :type alpha: float
    :param alpha: modifies the weight function W by changing the tradeoff between coverage and coverage overlap (default: 1.0).
    :type delta: float
    :param delta: number of delta allowed between any pair of gene sets (default: 0).
    :type lmbda: int
    :param lmbda: number of gene sets any gene can be a member of (default: 1).

    :returns: *p*-value (float) equal to the number of permuted matrices with greater *W'*, divided by the number of permuted matrices.

    **See also:** :func:`multi_dendrix.ILP`, :func:`load_permuted_matrices`.
    """
    count = 0.
    for mutation_data in permuted_matrices:
        multi_args = [ mutation_data, t, k_min, k_max, alpha, delta, lmbda ]
        collection_w_weights = Multi.ILP(*multi_args)
        if sum([w for gene_set, w in collection_w_weights]) >= W_prime:
            count += 1.
    return count / float(len(permuted_matrices))

def run(args):
    """Perform the matrix permutation test on the given collection of gene sets and the given directory of permuted mutation data."""
    # Load input
    if args.verbose: print 'Loading matrices...'
    W_prime = load_w_prime(args.collection_file)
    permuted_matrices = load_permuted_matrices(args.permuted_matrices_dir)

    # Evaluate the pathway set
    if args.verbose: print 'Evaluating input pathway set...'
    pval = matrix_permutation_test(W_prime, permuted_matrices, args.num_gene_sets,
		                           args.min_gene_set_size, args.max_gene_set_size,
		                           args.alpha, args.delta, args.lmbda)

    # Output and return results
    test_name = 'Matrix permutation pval: '
    if args.output_file:
        open(args.output_file, 'w').write(test_name + str(pval))
	
    if args.verbose: print test_name + str(pval)

if __name__ == "__main__": run(parse_args())