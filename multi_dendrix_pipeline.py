#!/usr/bin/python

"""
The Multi-Dendrix pipeline consists of the following steps:

1. Runs Multi-Dendrix on input mutation data for a range of number *t* of gene
   sets and maximum gene set size *kmax*. This yields a set of *collections* of
   gene sets (see :func:`multi_dendrix_pipeline.batch_multi_dendrix`). The
   motivation for running Multi-Dendrix on a range of parameter settings is that 
   while we expect genes in the same functional pathway to have approximate 
   exclusivity in their mutations, we do not know the sizes or the number of 
   these pathways. By running on a range of parameter sizes, we hope to 
   identify these true values (detailed in the next step).
2. Analyze the collections identified in 1. for the gene sets that appear 
   together consistently across parameter choices
   (see :func:`core_modules.extract_core_modules`). 
3. Evaluate each collection for statistical signficance using a *matrix 
   permutation test*.
   See :func:`multi_dendrix_pipeline.run_matrix_permutation_test` and
   :func:`matrix_permutation_test.matrix_permutation_test` for details.
4. Evaluate each collection for enrichment of protein-protein interactions on
   a protein-protein interaction (PPI) network. See
   :func:`multi_dendrix_pipeline.run_network_permutation_test` and
   :func:`network_tests.direct_interactions_test` for details.
5. Analyze all genes (or mutation classes) for (sub)type-specific mutations
   (if subtypes are known). See :func:`subtype_specific_genes.subtype_analysis`
   for details.
6. Output the results of the analysis as both text and HTML files.

At this time, the pipeline does not include any functions for preprocessing mutation data.
"""

# Load required modules and add the lib to the path
import sys, os
sys.path.insert(1, os.path.abspath('./lib'))
import multi_dendrix as Multi

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Runs Multi-Dendrix for a set of parameters. Evaluates the '\
                  'results and outputs them as text and as a website.'
    parser = argparse.ArgumentParser(description=description)

    # General options
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Name of output directory.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Flag verbose mode.')

    # Options for Multi-Dendrix
    parser.add_argument('-k_min', '--min_gene_set_size', required=True, type=int,
                        help='Minimum gene set size.')
    parser.add_argument('-k_max', '--max_gene_set_size', required=True, type=int,
                        help='Maximum gene set size.')
    parser.add_argument('-t_min', '--min_num_gene_sets', required=True, type=int,
                        help='Minimum number of gene sets.')
    parser.add_argument('-t_max', '--max_num_gene_sets', required=True, type=int,
                        help='Maximum number of gene sets.')
    parser.add_argument('-n', '--db_name', required=True,
                        help='Name of mutation data for use in output.')
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
    parser.add_argument('-a', '--alpha', type=float, default=1.0,
                        help='Parameter that changes weight function W by '\
                        'weighting the penalty of coverage overlap.')    
    parser.add_argument('--delta', type=int, default=0,
                        help='Number of overlaps allowed per gene set.')    
    parser.add_argument('--lmbda', type=int, default=1,
                        help='Number of gene sets a gene can be a member of.')
    
    # Options for core modules
    parser.add_argument('--stability_threshold', type=int, default=2,
                        help='Minimum number of gene sets two genes must '\
                             'both be a member of to be connected in the '\
                             'core modules.')

    # Options for (sub)type analysis
    parser.add_argument('--subtypes', default=False, action='store_true',
                        help='Perform (sub)type analysis.')
    parser.add_argument('--subtype_sig_threshold', default=0.05, type=float,
                        help='Significance threshold for subtype association '\
                             '(use 1.0 to report all associations).')

    # Options for permutation tests
    parser.add_argument('--network_test', default=False, action='store_true',
                        help='Perform network permutation test.')
    parser.add_argument('--weight_test', default=False, action='store_true',
                        help='Perform weight permutation test.')
    parser.add_argument('-ppi', '--network_edgelist', default=None,
                        help='PPI edgelist location.')
    parser.add_argument('--num_permuted_networks', default=5, type=int,
                        help='The number of permuted networks to create '\
                             '(only if a directory of permuted networks '\
                             'is not provided).')
    parser.add_argument('--permuted_networks_dir', default=None,
                        help='Directory of permuted networks.')
    parser.add_argument('--distance', default=False, action='store_true',
                        help='Flag average pairwise distance test.')
    parser.add_argument('--Q', default=100, type=int,
                        help='Multiplier of edge swaps for permuting networks.')

    parser.add_argument('--permuted_matrices_dir', default=None,
                        help='Directory of permuted matrices.')
    parser.add_argument('--num_permuted_matrices', default=5, type=int,
                        help='The number of permuted matrices to create '\
                             '(only if a directory of permuted matrices '\
                             'is not provided).')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args
 
def batch_multi_dendrix(args):
    """Runs Multi-Dendrix for each parameter setting on the input
    mutation data.

    **Returns:**
      A tuple containing the following:
        * **collections** (*dictionary*) - mapping of t -> k -> output of Multi-Dendrix
        * **mutation_data** (*tuple*) - mutation data tuple (see :func:`multi_dendrix.multi_dendrix` for details).
        * **runtime** (*float*) - total runtime (in seconds) of Multi-Dendrix on all the parameter settings
    """
    # Import required modules
    from time import time

    # Load mutation data used in each run
    start = time()
    include = Multi.white_and_blacklisting(args.patient_whitelist,
              args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, sample2include = include

    mutation_data = Multi.load_mutation_data_w_cutoff(args.mutation_matrix,
                    sample2include, gene2include, args.cutoff)
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    if args.verbose: print "* Mutation data: %s genes x %s patients" % (m, n)

    # Run Multi-Dendrix for the range of parameters 
    ts = range(args.min_num_gene_sets, args.max_num_gene_sets + 1)
    ks = range(args.min_gene_set_size, args.max_gene_set_size + 1)
    collections = dict( [(t, {}) for t in ts] )
    for t, k_max in [(t, k) for t in ts for k in ks]:
        multi_params = [ mutation_data, t, args.min_gene_set_size, k_max,
                         args.alpha, args.delta, args.lmbda ]
        collection_w_weights = Multi.ILP(*multi_params, verbose=args.verbose)
        collections[t][k_max] = zip(*collection_w_weights)

    return collections, mutation_data, time() - start

def run_network_permutation_test(args, collections, core_modules):
    """Runs the direct interactions or average pairwise distance test
    on each of the collections and the core_modules.

    **Returns**:
      * **evaluation** (*dictionary*) - a mapping of t -> k -> the network evaluation tuple of each collection (see :func:`network_tests.evaluate_collection` for details)
    """
    load_network = Multi.permute.ppi_network.load_network
    permute_network = Multi.permute.ppi_network.permute_network

    # Load original network and generate permuted networks
    G = load_network(args.network_edgelist)
    if args.permuted_networks_dir:
        network_files = [ args.permuted_networks_dir + "/" + fh
                          for fh in os.listdir(args.permuted_networks_dir)]
        Hs = [ load_network(H) for H in network_files ]
    else:
        Hs = [ permute_network(G, args.Q) for i in range(args.num_permuted_networks) ]

    # Perform network test
    evaluate_collection = Multi.evaluate.network.evaluate_collection

    evaluation = dict( [ (t, {}) for t in collections.keys() ] )
    for t in collections.keys():
        for k_max in collections[t].keys():
            gene_sets, weights = collections[t][k_max]
            results = evaluate_collection(gene_sets, G, Hs, args.distance)
            # test_name, statistic, pval, gene_set_results = results
            evaluation[t][k_max] = results

    evaluation["core_modules"] = evaluate_collection(core_modules, G, Hs,
                                                     args.distance)

    return evaluation

def run_matrix_permutation_test(args, collections, mutation_data):
    """Runs the matrix permutation test on each of the collections and the
       core_modules.

    **Returns**:
      * **evaluation** (*dictionary*) - a mapping of t -> k -> the network evaluation tuple of each collection (see :func:`network_tests.evaluate_collection` for details)
    """
    # Create shorter handles for important functions
    load_permuted_matrices = Multi.evaluate.matrix.load_permuted_matrices
    matrix_permutation_test = Multi.evaluate.matrix.matrix_permutation_test
    Permut = Multi.permute.mutation_data

    # Load / generate networks
    if args.permuted_matrices_dir:
        permuted_matrices = load_permuted_matrices(args.permuted_matrices_dir)
    else:
        m, n, genes, patients = mutation_data[:4]
        mutation2patients, patient2mutations = mutation_data[-2:]
        G = Permut.construct_mutation_graph(mutation2patients, patient2mutations)
        permuted_matrices = [ Permut.permute_mutation_data(G, genes, patients)
                              for i in range(args.num_permuted_matrices) ]

    # Perform network test
    evaluation = dict( [ (t, {}) for t in collections.keys() ] )
    for t in collections.keys():
        for k_max in collections[t].keys():
            gene_sets, weights = collections[t][k_max]
            test_args = [sum(weights), permuted_matrices,
                         t, args.min_gene_set_size, k_max, args.alpha,
                         args.delta, args.lmbda]
            pval = matrix_permutation_test(*test_args)
            evaluation[t][k_max] = pval

    return evaluation

def flatten_collections(collections):
    """Takes a dictionary of parameter settings to Multi-Dendrix results (as output by 
    :func:`batch_multi_dendrix`), and flattens the map into a list of collections."""
    all_collections = [ ]
    for t in collections.keys():
        for k_max in collections[t].keys():
            collection, weights = collections[t][k_max]
            all_collections.append( collection )
    return all_collections

def run(args):
    """Runs the whole :doc:`/pipeline` for the given command-line arguments."""
    # Run Multi-Dendrix for all parameter settings 
    collections, mutation_data, runtime = batch_multi_dendrix(args)
    
    # Extract the stable modules
    all_collections = flatten_collections(collections)
    core = Multi.core_modules.extract(all_collections, args.stability_threshold)
    core_modules, module_graph = core

    # Perform the permutation tests (if required)
    if args.weight_test:
        matrix_results = run_matrix_permutation_test(args, collections,
                                                     mutation_data)
    else: matrix_results = None

    if args.network_test:
        network_results = run_network_permutation_test(args, collections, 
                                                       core_modules)
    else: network_results = None

    evaluation  = network_results, matrix_results

    # Perform subtype analysis (if required)
    if args.subtypes and args.patient_whitelist:
        subtype_analysis = Multi.subtypes.subtype_analysis
        create_subtype_tbl = Multi.subtypes.create_subtype_tbl
        gene2specificity = subtype_analysis(mutation_data,
                                            args.patient_whitelist,
                                            args.subtype_sig_threshold)
        subtype_tbl = create_subtype_tbl(gene2specificity)

    else:
        gene2specificity, subtype_tbl = None, None
        if args.subtypes:
            print 'No patient whitelist (w/ (sub)types) provided, '\
                  'skipping (sub)type analysis.'

    # Create tables used for text and/or html output
    create_collection_tbls = Multi.output.create_collection_tbls
    create_params_tbl = Multi.output.create_params_tbl
    create_matrix_results_tbl = Multi.output.create_matrix_results_tbl
    create_network_results_tbl = Multi.output.create_network_results_tbl
    create_matrix_results_tbl = Multi.output.create_matrix_results_tbl

    collection_tbls = create_collection_tbls(args, collections,
                                               core_modules, evaluation)
    params_tbl = create_params_tbl(args, mutation_data)

    if args.network_test: network_tbl =  create_network_results_tbl(evaluation[0])
    else: network_tbl = None
    if args.weight_test: matrix_tbl = create_matrix_results_tbl(evaluation[1])
    else: matrix_tbl = None

    # output results to text and html
    text_output_args = [ args, collection_tbls, runtime, params_tbl,
                         network_tbl, matrix_tbl, subtype_tbl ]
    Multi.output.output_to_text(*text_output_args)

    Multi.output.output_to_html(args, collections, runtime, module_graph,
        evaluation, params_tbl, subtype_tbl, gene2specificity)

if __name__ == "__main__": run(parse_args())
