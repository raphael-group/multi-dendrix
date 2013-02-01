#!/usr/bin/env python

# Import required modules
import sys
try: import cplex
except ImportError:
    print "Couldn't import CPLEX. Check your environment's PYTHONPATH "\
          "variable. For details on the CPLEX python module, please visit "\
          "http://bit.ly/KL7PVc."
    sys.exit(1)

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Runs Multi-Dendrix to find the optimal set of t gene sets '\
                  'of k genes for the weight function W\'.'
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
    parser.add_argument('-t', '--num_gene_sets', required=True, type=int,
                        help='Desired number of gene sets.')
    parser.add_argument('-k_min', '--min_gene_set_size', type=int,
                        help='Minimum gene set size.')
    parser.add_argument('-k_max', '--max_gene_set_size', type=int,
                        help='Maximum gene set size.')
    parser.add_argument('-k', '--gene_set_size', type=int, default=None,
                        help='Gene set size.')
    parser.add_argument('-a', '--alpha', type=float, default=1.0,
                        help='Parameter that changes weight function W by chang'\
                        'ing the tradeoff between coverage and exclusivity')
    parser.add_argument('-o', '--output_file', default=None,
                        help='Name of output file.')    
    parser.add_argument('-d', '--delta', type=int, default=0,
                        help='Number of overlaps allowed per gene set.')    
    parser.add_argument('-l', '--lmbda', type=int, default=1,
                        help='Number of gene sets a gene is allowed to be in.') 
    parser.add_argument('--max_time', type=int, default=1e75,
                        help='Maximum amount of time to be spent on an optimization.')   
    parser.add_argument('--max_mem', type=int, default=1e75,
                        help='Maximum tree memory to be used by CPLEX.')
    parser.add_argument('-q', '--quiet', default=False, action="store_true",
                        help='Quiet output flag.')
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    return args


########################## Multi-Dendrix ILP ###################################
# Helper functions for the multi_dendrix ILP function
def A(patient2mutations, gene, patient):
    '''Mutation  matrix accessor function (see :func:`ILP`).

    :type patient2mutations: dictionary
    :param patient2mutations: mapping of patients to the genes in which they have mutations.
    :type gene: string
    :param gene: gene name.
    :type patient: string
    :param patient: patient name.
    :returns: 1 if gene g is mutated in the given patient, 0 otherwise.

    **Examples:**
      >>> A({"TCGA-01" : ["G2", "G3"], "G2", "TCGA-01"})
      1
      >>> A({"TCGA-01" : ["G2", "G3"], "G1", "TCGA-01"})
      0
    '''
    return int(gene in patient2mutations[patient])

def var_to_gene(var):
    '''Extracts the gene name from the variable name used in the ILP'''
    if "->" in var: return var.split("->")[0]
    else: return var
    
def path_num(var):
    '''Extracts the set number from a variable name, and corrects the index'''
    if "->" in var:
        try: return int(var.split("->")[1]) - 1
        except ValueError: return 0
    else: return 0

# CPLEX memory options (see http://goo.gl/dJV6F)
MEMORY_ONLY        = 1
ON_DISK_COMPRESSED = 3

# Main function for Multi-Dendrix ILP
def ILP(mutation_data, t=2, k_min=2,
                  k_max=3, alpha=1.0, delta=0.0, lmbda=1.0, verbose=False,
                  max_mem=None, max_time=None):
    '''Implementation of Multi-Dendrix ILP. Sets up ILP, uses CPLEX Python to solve it, and parses the results.

    :type t: int
    :param t: number of gene sets (default: 2).
    :type k_min: int
    :param k_min: minimum gene set size (default: 2).
    :type k_max: int
    :param k_max: maximum gene set size (default: 3).
    :type alpha: float
    :param alpha: modifies the weight function W by changing the tradeoff between coverage and coverage overlap (default: 1.0).
    :type delta: float
    :param delta: number of delta allowed between any pair of gene sets (default: 0).
    :type lmbda: int
    :param lmbda: number of gene sets any gene can be a member of (default: 1).
    :param verbose: outputs progress to stdout (default: False).
    :type max_mem: float
    :param max_mem: amount of memory to constrain the CPLEX optimizer to (default: unconstrained).
    :type max_time: float
    :param max_time: amount of time to constrain the CPLEX optimizer to (default: unconstrained).

    :returns: A list of *t* tuples sorted by weight W, where each tuple contains 1) a gene set that is part of an optimal solution for the given input data, and 2) the weight *W* of that gene set (recall that the function *W* changes with the parameter *alpha*).
    
    **Examples:**
      A view of example input:
        >>> mutation2patients = {"G1" : ["TCGA-01", "TCGA-03"], "G2" : ["TCGA-02"],
        "G3" : ["TCGA-01", "TCGA-02", "TCGA-03"], "G4" : ["TCGA-02", "TCGA-04"]}
        >>>> patient2mutations = {"TCGA-01" : ["G1", "G3"], "TCGA-02" : ["G2", "G2"],
                                  "TCGA-03" : ["G1", "G3"], "TCGA-04" : ["G4"]}
        >>> genes = ["G1", "G2", "G3", "G4"]
        >>> patients = ["TCGA-01", "TCGA-02", "TCGA-03", "TCGA-04"]
        >>> data = (4, 4, genes, patients, mutation2patients, patient2mutations)
      The results of Multi-Dendrix run for *t* = *kmax* = *kmin* = 2:
        >>> ILP(data, 2, 2, 2)
        [(['G3', 'G4'], 3), (['G1', 'G2'], 3)]

    **See also:** :func:`load_mutation_data`.
    '''
    # Parse mutation data into shorter variable handles
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

    # Default parameters for ILP
    try: ilp, my_obj, my_sense, my_rhs = cplex.Cplex(), [], "L", [1]
    except cplex.exceptions.CplexSolverError:
        print 'Exiting because of CPLEX license issue.'
        sys.exit(2)

    # Set params for ILP
    ilp.set_results_stream(None) # Change this to have CPLEX output to stdout
    ilp.parameters.mip.strategy.file.set(MEMORY_ONLY)
    #if max_time:
    #    ilp.parameters.timelimit.set(max_time)
    #    if verbose: print 'Time limit set to', max_time
    if max_mem:
        ilp.parameters.mip.limits.treememory.set(max_mem)
        if verbose: print 'Max tree memory set to', max_mem

    patients = patient2mutations.keys()

    # The ILP needs a variable for each gene in each gene set it could be a
    # member of. Creates variable names of the form  {gene}-{set number}
    # |all_names| = |genes| x |gene_sets|
    l = map(lambda g: [g + "->" + str(i) for i in range(1, t + 1)], genes)
    all_names, over_is = [item for sublist in l for item in sublist], []

    # initialize r variables
    rs = ["rvariable" + str(i) + "-" + str(j) for i in range(len(patients))
          for j in range(1, t + 1)]
    rs = []
    for j in range(1, t + 1):
        rs.append(["rvariable-" + str(i) + "-" + str(j)
                   for i in range(len(patients))])
    flat_rs =  [item for sublist in rs for item in sublist]
    rs_coefs = [1.0 + alpha] * len(flat_rs)
    C = 65536

    # Encode W' as the objective function
    for v in all_names:
        my_obj.append(len(mutation2patients[var_to_gene(v)]) * -alpha)

    # append r vars to the set of variables and r vars' coefficients to the obj
    # function. genes' vars are already initialized to |mutation2patients[g]|
    my_obj += rs_coefs
    all_names += flat_rs

    if t > 1:
        for i in range(1, t + 1):
            for j in range(i + 1, t + 1):
                over_is += [g + ".." + str(i) + ".." + str(j) for g in genes]
        all_names += over_is
        my_obj += [0.0] * len(over_is)

    # initialize vars and their bounds {0,1}
    lbs, ubs = [0] * len(all_names), [1] * len(all_names)
    my_types = [ilp.variables.type.integer] * len(all_names)
    ilp.variables.add(names = all_names, ub = ubs, lb = lbs,
                      types = my_types, obj = my_obj)
    ilp.objective.set_sense(ilp.objective.sense.maximize)

    # add gene set size constraints, if any
    for i in range(1, t + 1):
        my_names = [g + "->" + str(i) for g in genes]
        size_coefs = [1.0] * len(my_names)
        if k_max != 0 and k_min == k_max:
            ilp.linear_constraints.add(senses = "E", rhs = [k_min],
                                       names = ["eq" + str(i)])
            ilp.linear_constraints.set_linear_components(
                [("eq" + str(i), [my_names, size_coefs])])
        elif type([]) == type(k_max):
            ilp.linear_constraints.add(senses = "E", rhs = [k_max[i-1]],
                                       names = ["eq" + str(i)])
            ilp.linear_constraints.set_linear_components(
                [("eq" + str(i), [my_names, size_coefs])])  
        else:
            ilp.linear_constraints.add(senses = "L", rhs = [k_max],
                                       names = ["max" + str(i)])
            ilp.linear_constraints.set_linear_components(
                [("max" + str(i), [my_names, size_coefs])])
            ilp.linear_constraints.add(senses = "G", rhs = [k_min],
                                       names = ["min" + str(i)])
            ilp.linear_constraints.set_linear_components(
                [("min" + str(i), [my_names, size_coefs])])

    # add non-overlapping gene sets constraint
    for g in genes:
        my_names = [g + "->"  + str(i) for i in range(1, t + 1)]
        coefficients = [1] * len(my_names)
        ilp.linear_constraints.add(senses = my_sense, rhs = [lmbda],
                                   names = [g])
        ilp.linear_constraints.set_linear_components([(g, [my_names,
                                                           coefficients])])

    # add r constraints
    for i in range(len(patients)):
        coefficients = [A(patient2mutations, g, patients[i]) for g in genes]
        my_coefficients = coefficients + [-1.0]
        for j in range(1, t + 1):
            r = rs[j-1][i]
            my_names = [g + "->" + str(j) for g in genes] + [r]
            ilp.linear_constraints.add(senses = "G", rhs = [0.0],
                                       names = [str(i) + str(j) + "-rn"])
            ilp.linear_constraints.set_linear_components(
                [(str(i) + str(j) + "-rn", [my_names, my_coefficients])])

    # allow at most one overlap per gene set, if there's more than one gene set
    ###########################################################################
    if t > 1 and delta > 0:
        # first add pairwise gene set indicators
        for ind in over_is:
            g, i, j = ind.split("..")[0], ind.split("..")[1], ind.split("..")[2]
            my_names = [g + "->"  + i, g + "->" + j, ind]
            coefficients = [1.0, 1.0, -1.0 * (C+1)]
            ilp.linear_constraints.add(senses = "G", rhs = [C * -1.0 + 1.0],
                                       names = [g + i + j + "1"])
            ilp.linear_constraints.set_linear_components(
                [(g + i + j + "1", [my_names, coefficients])])
            coefficients = [1.0, 1.0, -1.0 * C]
            ilp.linear_constraints.add(senses = "L", rhs = [1.0],
                                       names = [g + i + j + "2"])
            ilp.linear_constraints.set_linear_components(
                [(g + i + j + "2", [my_names, coefficients])])

        # then constrain pairwise gene set indicators so only one gene can
        # be shared between two gene sets
        for i in range(1, t + 1):
            for j in range(i + 1, t + 1):
                my_names = [g + ".." + str(i) + ".." + str(j) for g in genes]
                coefficients = [1.0] * len(my_names)
                ilp.linear_constraints.add(senses = "L", rhs = [delta],
                                           names = ["k" + str(i) + str(j)])
                ilp.linear_constraints.set_linear_components(
                    [("k" + str(i) + str(j), [my_names, coefficients])])
                
    ###########################################################################

    # Solve ILP and convert solution back to genes
    ilp.solve()
    ans = ilp.solution.get_values()
    collection = [[] for i in range(t)]
    for var, val in zip(all_names, ans):
        if str(val) == "1.0" or float(val) > 0.5:
            # be sure to ignore dummy variables
            if var.count("rvariable") == 0 and var.count("..") == 0:
                collection[path_num(var)].append(var_to_gene(var))
    
    # sort each gene set by the highest to lowest coverage genes
    for p in collection: p.sort(key=lambda g: len(mutation2patients[g]),
                                reverse=True)
    
    # sort pathways by weight (highest first) for output
    collection_with_weights = [(p, W(mutation2patients, p, alpha))
                               for p in collection]
    collection_with_weights.sort(key=lambda (p, w): w)
    collection_with_weights.reverse()

    # Output run summary if verbose
    if verbose:
        print "done!\n"
        print "------------- Multi-Dendrix Run Summary ----------------"
        print "Num gene sets:", t, "\nMin:", k_min, "\nMax:", k_max
        print "delta:", delta, "\nMax Sets / Gene:", lmbda
        print "Alpha:", alpha
        print "\nCollection:"
        W_prime = 0
        for i in range(1, len(collection_with_weights) + 1):
            p, w = collection_with_weights[i-1]
            print "\t", ", ".join(p)
            W_prime += w
            print "\tW(M) =", w
        print "W'(M):", W_prime
        print "\n-----------------------------------------------------\n"

    return collection_with_weights

################################################################################

####### Helper functions for running Multi-Dendrix and analyzing results #######
# Weight function from Dendrix (Vandin et al. (2012) Genome Research 22:375-85).
def W(mutation2patients, gene_set, alpha):
    """Calculates the weight W of a set of genes, defined as the weighted
    difference between their total coverage and coverage overlap.

    :type mutation2patients: dictionary
    :param mutation2patients: mapping gene / mutation class to the patients with mutations in that gene.
    :type gene_set: list
    :param gene_set: a list of gene names.
    :type alpha: float
    :param alpha: Multi-Dendrix parameter alpha (see :func:`ilp` for reference).
    :returns: weight *W*

    **Example:**

    >>> h = { "g1" : [3, 4, 5], "g2" : [1, 2, 3], "g3" : [5, 6, 7]}
    >>> gset = ["g1", "g2", "g3"]
    >>> W(h, gset, 2)
    5
    """
    patients, coverage_overlap = set(), 0
    for gene in gene_set:
        targets  = mutation2patients[ gene ]
        patients.update( targets )
        coverage_overlap += len( targets )

    coverage = len(patients)
    return int( (1+alpha) * coverage - alpha * coverage_overlap )

def load_mutation_data_w_cutoff(file_loc, patient_wlst=None,
                                gene_wlist=None, cutoff=0):
    """Loads the mutation data in the given file, restricting to genes with a given mutation frequency. 

    :type file_loc: string
    :param file_loc: Location of mutation data file.
    :type patient_wlst: dictionary
    :param patient_wlst: Maps patient IDs to whether they should be included in the analyzed mutation data.
    :type gene_wlist: dictionary
    :param gene_wlist: Maps genes to whether they should be included in the analyzed mutation data.
    :type cutoff: int
    :param cutoff: Minimum mutation frequency a gene must have to be included in the analyzed mutation data.
    :returns: Mutation data tuple (see :func:`load_mutation_data`).

    **Example:**
      A view into the example data:
        >>> file_loc = 'test.m2'
        >>> open(file_loc).readlines() # view of the  data
        ["TCGA-01\\tG1\\tG3\\tG5\\n", "TCGA-02\\tG2\\tG1\\tG4\\n"]
    
      Mutation data with no cutoff:
        >>> load_mutation_data_w_cutoff(file_loc, cutoff=0) # mutation data with no cutoff
        (5, 2, ["G1", "G2", "G3", "G4", "G5"], ["TCGA-01", "TCGA-02"],
            {"G1" : ["TCGA-01", "TCGA-02"], "G2" : ["TCGA-02"], "G3" : ["TCGA-01"],
            "G4" : ["TCGA-2"], "G5" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G3", "G5"], "TCGA-02" : ["G2", "G1", "G4"]})
      
      Mutation data with a cutoff of 2:
        >>> load_mutation_data_w_cutoff(file_loc, cutoff=2)
        (1, 2, ["G1"], ["TCGA-01", "TCGA-02"], {"G1" : ["TCGA-01", "TCGA-02"]},
            {"TCGA-01" : ["G1"], "TCGA-02" : ["G1"]})
    
    **See also:**
    :func:`white_and_blacklisting`, :func:`load_mutation_data`

    """
    mutation_data = load_mutation_data(file_loc, patient_wlst, gene_wlist)
    if cutoff == 0: return mutation_data
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    for gene in genes:
        if len(mutation2patients[gene]) < cutoff:
            for patient in mutation2patients[gene]:
                patient2mutations[patient].remove( gene )
            del mutation2patients[gene]
    
    m, genes = len(mutation2patients.keys()), mutation2patients.keys()
    return m, n, genes, patients, mutation2patients, patient2mutations

def load_mutation_data(file_loc, patient_whitelist=None, gene_whitelist=None):
    """Loads the mutation data in the given file. 

    :type file_loc: string
    :param file_loc: Location of mutation data file.
    :type patient_whitelist: dictionary
    :param patient_whitelist: Maps patient IDs to whether they should be included in the analyzed mutation data.
    :type gene_whitelist: dictionary
    :param gene_whitelist: Maps genes to whether they should be included in the analyzed mutation data.
    :rtype: Tuple

    **Returns:**
      * **m** (*int*) - number of patients.
      * **n** (*int*) - number of genes.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **mutation2patients** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patient2mutations** (*dictionary*) - mapping of patients to the genes they have mutated.

    **Example:**
      A view into example data:
        >>> file_loc = 'test.m2'
        >>> open(file_loc).readlines()
        ["TCGA-01\\tG1\\tG3\\tG5\\n", "TCGA-02\\tG2\\tG1\\tG4\\n"]
      
      Mutation data with no whitelisting:
        >>> load_mutation_data(file_loc)
        (5, 2, ["G1", "G2", "G3", "G4", "G5"], ["TCGA-01", "TCGA-02"],
            {"G1" : ["TCGA-01", "TCGA-02"], "G2" : ["TCGA-02"], "G3" : ["TCGA-01"],
            "G4" : ["TCGA-2"], "G5" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G3", "G5"], "TCGA-02" : ["G2", "G1", "G4"]})
      
      Mutation data with patient whitelisting only:
        >>> patient_wlst = {"TCGA-01" : False, "TCGA-02" : True}
        >>> load_mutation_data_w_cutoff(file_loc, patient_wlst)
        (3 , 1, ["G1", "G2", "G4"], ["TCGA-02"],
            {"G1" : ["TCGA-02"], "G2" : ["TCGA-02"], "G4" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G2", "G4"]})
      
      Mutation data with patient and gene whitelisting:
        >>> gene_wlst = {"G1" : True, "G2" : True, "G4" : False}
        >>> load_mutation_data_w_cutoff(file_loc, patient_wlst, gene_wlst
        (2 , 1, ["G1", "G2"], ["TCGA-02"], {"G1" : ["TCGA-02"], "G2" : ["TCGA-01"]},
            {"TCGA-01" : ["G1", "G2"]})

    **See also:**
    :func:`white_and_blacklisting`, :func:`load_mutation_data_w_cutoff`
    """
    patient2mutations, mutation2patients = {}, {}
    for line in [l.rstrip() for l in open(file_loc) if not l.startswith('#')]:
        arr = filter(lambda g: g != "", line.split('\t'))
        patient, mutations = arr[0], arr[1:]

        if gene_whitelist:
            mutations = filter(lambda g: gene_whitelist[g], mutations)
        if patient_whitelist and not patient_whitelist[patient]: continue
        
        patient2mutations[patient] = mutations
        for gene in mutations:
            if not gene in mutation2patients.keys():
                mutation2patients[gene] = set([patient])
            else:
                mutation2patients[gene].add(patient)

    genes, patients = mutation2patients.keys(), patient2mutations.keys()
    m, n = len(genes), len(patients)
    return m, n, genes, patients, mutation2patients, patient2mutations

def white_and_blacklisting(patient_wlst=None, patient_blst=None, gene_wlst=None,
                           gene_blst=None):
    '''Reconciles the different white- and blacklists provided as input into Multi-Dendrix.

    :type patient_wlst: string
    :type patient_blst: string
    :type gene_wlst: string
    :type gene_blst: string
    :param patient_wlst: File location of patients to be *included* in analyzed mutation data.
    :param patient_blst: File location of patients to be *excluded* in analyzed mutation data.
    :param gene_wlst: File location of patients to be *included* in analyzed mutation data.
    :param gene_blst: File location of patients to be *excluded* in analyzed mutation data.

    **Returns:**
      * gene2include (*dictionary*): mapping of genes to whether they should be included in the analyzed mutation data.
      * patient2include (*dictionary*): mapping of patients to whether they should be included in the analyzed mutation data.

    **Examples:**
      *(For brevity, examples are for patient white- and blacklists only)*

      Patient whitelisting only:
        >>> patient_wlst = 'patient.wlst'
        >>> open(patient_wlst).readlines()
        ["TCGA-01", "TCGA-02", "TCGA-03"]
        >>> white_and_blacklisting(patient_wlst)
        (defaultdict(<function <lambda>>, {}), {"TCGA-01", "TCGA-02", "TCGA-03"})
      Conflicting patient white- and blacklists (whitelist has final word):
        >>> patient_blst = 'patient.blst'
        >>> open(patient_wlst).readlines()
        ["TCGA-02", "TCGA-04"]
        >>> white_and_blacklisting(patient_wlst)
        (defaultdict(<function <lambda>>, {}), {"TCGA-01", "TCGA-02", "TCGA-03"})


    **See also:** :func:`load_mutation_data`, :func:`load_mutation_data_w_cutoff`.

    '''

    # Blacklisting and whitelisting works as follows. If a whitelist is passed in,
    # only genes/patients in the whitelist and NOT in the blacklist are allowed. If
    # only a blacklist is passed in, all genes/patients not on the blacklist are
    # allowed.
    from collections import defaultdict
    if patient_wlst:
        patient2include = defaultdict(lambda : False)
        patient_whitelist = [l.rstrip().split()[0] for l in open(patient_wlst)
                             if not l.startswith("#")]

        for p in patient_whitelist: patient2include[p] = True
    else:
        patient_whitelist = None
        patient2include = defaultdict(lambda : True)

    if patient_blst:
        patient_blacklist = [l.rstrip() for l in open(patient_blst)]
        for p in patient_blacklist: patient2include[p] = False

    if gene_wlst:
        gene2include = defaultdict(lambda : False)
        gene_whitelist = set([l.split()[0] for l in open(gene_wlst)])
        for g in gene_whitelist: gene2include[g] = True
    else: gene_whitelist, gene2include = None, defaultdict(lambda : True)
    
    if gene_blst:
        gene_blacklist = [l.rstrip() for l in open(gene_blst)]
        for g in gene_blacklist: gene2include[g] = False

    return gene2include, patient2include

def run(args):
    """Runs Multi-Dendrix for the given arguments, and outputs the results to file."""
    # Reconcile the input white- and blacklists
    include = white_and_blacklisting(args.patient_whitelist,
              args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, patient2include = include

    mutation_data = load_mutation_data_w_cutoff(args.mutation_matrix,
                                               patient2include, gene2include,
                                               args.cutoff)
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

    # Run Multi-Dendrix
    if args.gene_set_size:
        k_min, k_max = args.gene_set_size, args.gene_set_size,
    else:
        k_min, k_max = args.min_gene_setsize, args.max_gene_set_size
    
    collection_with_weights = ILP(mutation_data, args.num_gene_sets,
        k_min, k_max, args.alpha, args.delta, args.lmbda, max_mem=args.max_mem,
        max_time=args.max_time, verbose=not args.quiet)

    # Output results
    rows = [str(weight) + ' ' + 
            ' '.join([g + '[' + ','.join(mutation2patients[g]) + ']'
                      for g in gene_set])
            for gene_set, weight in collection_with_weights]
    if args.output_file:
        open(args.output_file, 'w').write('\n'.join(rows))
    return collection_with_weights

if __name__ == "__main__": run(parse_args())
