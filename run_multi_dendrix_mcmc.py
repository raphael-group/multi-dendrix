#!/usr/bin/python

from multi_dendrix_mcmc import *
import multi_dendrix as Multi

# Create shorter handles for important functions
load_permuted_matrices = Multi.evaluate.matrix.load_permuted_matrices
matrix_permutation_test = Multi.evaluate.matrix.matrix_permutation_test
Permut = Multi.permute.mutation_data

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Runs Multi-Dendrix MCMC.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-c', '--cutoff', type=int, default=0, 
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-g', '--gene_whitelist', default=None,
                        help='File of genes to be included.')
    parser.add_argument('-bg', '--gene_blacklist', default=None,
                        help='File of genes to not be included.')
    parser.add_argument('-ws', '--sample_whitelist', default=None,
                        help='File of samples to be included.')
    parser.add_argument('-bs', '--sample_blacklist', default=None,
                        help='File of samples to not be included.')
    parser.add_argument('-k', '--gene_set_size', type=int, required=True,
                        help='Gene set size.')
    parser.add_argument('-t', '--num_gene_sets', type=int, required=True,
                        help='Number of gene sets.')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Path to output file.')    
    parser.add_argument('-q', '--quiet', default=False, action="store_true",
                        help='Quiet output flag.')
    parser.add_argument('-n', '--num_iterations', type=int, default=pow(10, 6),
                        help='Number of iterations of MCMC.')
    parser.add_argument('-s', '--step_length', type=int, default=100,
                        help='Number of iterations between samples.')
    parser.add_argument('-p', '--num_permutations', type=int, default=0,
                        help='Number of permutations to perform.')
    parser.add_argument('--permuted_matrices_dir', default=None,
                        help='Directory of permuted matrices.')
    parser.add_argument('--parallel', action='store_true', default=False,
                        help='Use all available cores for permutation test '\
                              '(recommended for machines with at least 8 cores.')
    
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)

    return args

###############################################################################
# Convert mutation data to integers for easy C manipulation

def convert_mutation_data((m, n, genes, samples, gene2cases, sample2genes)):
    # Create a gene2index and index2gene mapping
    genes, samples = list(genes), list(samples)
    gene2index = dict([ (genes[i], i) for i in range(m) ])
    index2gene = dict([ (i, genes[i]) for i in range(m) ])

    # Create a dictionary of gene indices to sample indices
    gene2case_index = [ [samples.index(s) for s in gene2cases[genes[i]]]
                        for i in range(m)]

    return gene2index, index2gene, gene2case_index


def weight(genes, gene2cases):
    coverage = [ s for g in genes for s in gene2cases[g] ]
    return 2*len(set(coverage)) - len(coverage)

def convert_solns(index2gene, solns):
    new_solns = []
    for arr in solns:
        arr.sort(key=lambda M: M[-1], reverse=True)
        S = tuple( frozenset([index2gene[g] for g in M[:-1] ]) for M in arr )
        W = [ int(M[-1]) for M in arr ]
        new_solns.append( (S, W) )
    return new_solns

def collate_results( solns ):
    FREQ = 0
    results = dict()
    for S, Ws in solns:
        if S in results: results[S]["freq"] += 1
        else: results[S] = {"Ws": Ws, "freq": 1, "W_prime": sum(Ws)}
    return results

def run_matrix_permutation_test(npermutations, permuted_matrices_dir, mutation_data,
                                k, t, N, s, flat_results, parallel):
    # Find the max weight gene set from the MCMC on permuted data
    if npermutations > 0: print "* Running permutation test..."
    else: return {}

    # Load / generate permuted matrices. Otherwise, create permutation graph
    if permuted_matrices_dir:
        permuted_matrices = load_permuted_matrices(permuted_matrices_dir)
        mutation_graph, genes, patients = None, None, None
    else:
        permuted_matrices = [ None ] * npermutations
        m, n, genes, patients = mutation_data[:4]
        gene2cases, sample2genes = mutation_data[-2:]
        mutation_graph = Permut.construct_mutation_graph(gene2cases, sample2genes)

    # 
    permute_pval = dict()
    if parallel:
        import multiprocessing as mp
        pool         = mp.Pool()
        permute_args = [ (mcmc, k, t, N, s, mutation_graph,
                          genes, patients, permuted_matrices[i])
                        for i in range(npermutations) ]
        max_weights  = pool.map( permutation_test_wrapper, permute_args )
        pool.close()
        pool.join()
    else:
        max_weights = [permutation_test( mcmc, k, t, N, s,
                       mutation_graph, genes, patients, permuted_matrices[i])
                       for i in range(npermutations)]

    # Compute the significance of each gene set
    for i, (S, freq, W_prime, Ws) in enumerate(flat_results):
        # Only evalaute the top 10 gene sets
        if i >= 25: permute_pval[S] = "--"
        
        # Else P-values are the number of permutations for which the max weight
        # gene set found had a weight less than the current gene set's weight
        else:
            more_extreme = len([Ws[i] for Ws in max_weights if Ws[i] > W_prime ])
            if more_extreme == 0:
                permute_pval[S] = "P<%s" % (1./npermutations)
            else:
                permute_pval[S] = "P=%s" % (float(more_extreme)/npermutations)

    return permute_pval

# Returns the top scoring set after permuting the given mutation data
def permutation_test( mcmc, k, t, N, s, mutation_graph, genes, patients, permuted_matrix):
    # Permute the mutation data (if necessary)
    if not permuted_matrix:
        permuted_matrix = Permut.permute_mutation_data(mutation_graph, genes, patients, Q=1)
    m, n, genes, samples, gene2cases, sample2genes = permuted_matrix

    # Convert the mutation data to integers for easy C parsing
    gene2index, index2gene, gene2case_index = convert_mutation_data(permuted_matrix)

    # Run the MCMC algorithm with the given parameters and extract the top weight
    solns = mcmc(k, t, m, n, gene2case_index, N, s, 0)
    solns_w_weights = convert_solns( index2gene, solns )    

    return sorted([ sum(arr[-1]) for arr in solns ])[:25]

# Wrapper for permutation test that only takes one argument
# (used for multiprocessing module)
def permutation_test_wrapper( args ): return permutation_test(*args)

def run( args ):
    # Parse args into shorter variable handles
    mutation_matrix = args.mutation_matrix
    cutoff, gene_whitelist = args.cutoff, args.gene_whitelist
    k, t = args.gene_set_size, args.num_gene_sets
    N, s = args.num_iterations, args.step_length
    npermutations = args.num_permutations

    # Load mutation data
    print "* Loading mutation data..."
    include = Multi.white_and_blacklisting(args.sample_whitelist,
              args.sample_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, sample2include = include

    mutation_data = Multi.load_mutation_data_w_cutoff(args.mutation_matrix,
                                               sample2include, gene2include,
                                               args.cutoff)
    m, n, genes, patients, gene2cases, sample2genes = mutation_data

    # Convert the mutation data to integers for easy C parsing
    gene2index, index2gene, gene2case_index = convert_mutation_data(mutation_data)

    # Run the MCMC algorithm
    print "* Running MCMC algorithm..."
    solns = mcmc(k, t, m, n, gene2case_index, N, s, int(not args.quiet))

    # Collate the results and sort them descending by sampling frequency
    solns_w_weights = convert_solns( index2gene, solns )
    results         = collate_results( solns_w_weights )

    # Flatten the results
    print "* Flattening the results..."
    flat_results = [ [S, data["freq"], data["W_prime"], data["Ws"]] for S, data in results.items() ]
    flat_results.sort(key=lambda (S, F, W_prime, Ws): F, reverse=True)

    # Perform the permutation test (if necessary)
    if args.num_permutations > 0:
        permute_pval = run_matrix_permutation_test(npermutations, args.permuted_matrices_dir,
                                        mutation_data, k, t, N, s, flat_results, args.parallel)

    # Output results
    print "* Outputting results..."
    if args.num_permutations > 0:
        output = ["#Sampling Frequency\tWeight W'\tP-value ({} permutations)"\
                  "\tGene set\tWeight W".format(args.num_permutations)]
    else:
        output = ["#Sampling Frequency\tWeight W'\tGene set\tWeight W"]
    
    for S, freq, W_prime, Ws in flat_results:
        for i, M in enumerate(S):
            if i == 0:
                row = [ freq, W_prime, ", ".join(sorted(M)), Ws[i] ]
                if args.num_permutations > 0:
                    row.insert(2, permute_pval[S])
            else:
                row = [ "", "", "", ", ".join(sorted(M)), Ws[i] ]            

            output.append( "\t".join(map(str, row)) )

    open(args.output_file, "w").write( "\n".join(output) )

if __name__ == "__main__": run( parse_args() )
