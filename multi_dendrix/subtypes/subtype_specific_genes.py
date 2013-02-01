#!/usr/bin/python

# Load required modules
# Try and load scipy's fisher's exact test function
try:
	from scipy.stats import fisher_exact as pvalue
	def fisher_exact(tbl): return pvalue(tbl, alternative='greater')[1]
except ImportError:
    try:
        from fisher import pvalue as pvalue
        def fisher_exact(tbl): return pvalue(*tbl).right_tail
    except ImportError:
        import sys
        print 'Fatal Error: Neither SciPyv0.11 or fisher0.1.4 modules '\
              '(http://goo.gl/zYrLr) are installed.'
        sys.exit(1)

def parse_args(input_list=None):
    # Parse arguments
    import argparse
    class Args: pass
    args = Args()
    description = 'Calculates whether any genes are subtype specific for the '\
                   'given mutation data.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-c', '--cutoff', type=int, default=0, 
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-p', '--patient_whitelist', required=True,
                        help='Space-separated file of patient IDs and their '\
                              '(sub)type to be tested against.')
    parser.add_argument('-bp', '--patient_blacklist', default=None,
                        help='File of patients to be excluded.')
    parser.add_argument('-g', '--gene_whitelist', default=None,
                        help='File of genes to be included.')
    parser.add_argument('-bg', '--gene_blacklist', default=None,
                        help='File of genes to be excluded.')
    parser.add_argument('-o', '--output_file', default=None,
                        help='Name of output file.')
    parser.add_argument('--sig_threshold', default=0.05, type=float,
                        help='Significance threshold.')
    parser.add_argument('-a', '--all', default=False, action='store_true',
                        help='Flag to output all associations.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Flag verbose mode.')

    # If called from the command line, parse command line args.
    if input_list: parser.parse_args(input_list, namespace=args)
    else: parser.parse_args(namespace=args)
    
    return args	

def ty_contingency_table(ty, ty2mutations, tys, ty2numpatients):
    """Constructs the contigency table used by Fisher's exact test for subtype-specific mutation analysis.
    
    Contingency table:
      .. image:: /_static/subtype_specific_genes_cont_table.*

    :type ty: string
    :param ty: Name of (sub)type.
    :type ty2mutations: dictionary
    :param ty2mutations: Mapping of each subtype to the mutations (in one particular gene) in that subtype.
    :type tys: list of strings
    :param tys: List of all (sub)types.
    :type ty2numpatients: dictionary
    :param ty2numpatients: mapping of each (sub)type to the number of patients of that subtype.

    :returns: The tuple (*x11*, *x10*, *x01*, *x00*).

    **Examples:**
      A view of the input data:
        >>> ty2mutations = {"Luminal A" : ["TCGA-01", "TCGA-02"], "Luminal B" : []}
        >>> ty2numpatients = {"Luminal A" : 2, "Luminal B" : 1}
        >>> tys = ["Luminal A", "Luminal B"]
      A simple example:
        >>> ty_contingency_table("Luminal A", ty2mutations, tys, ty2numpatients)
        (2, 0, 0, 1)
    **See also:** :func:`subtype_specificity`, :func:`subtype_analysis`.
    """
    type_mutations = len(ty2mutations[ty])
    non_type_mutations = sum([len(ty2mutations[ty2]) for ty2 in tys if ty != ty2])

    type_normal = ty2numpatients[ty] - len(ty2mutations[ty])
    non_type_normal = sum([ty2numpatients[ty2] - len(ty2mutations[ty2]) for ty2 in tys if ty != ty2])

    return type_mutations, type_normal, non_type_mutations, non_type_normal

def subtype_specificity(gene, patient2ty, ty2numpatients, mutation2patients):
    """Performs a statistical test on a gene for each given (sub)type.

    **Test**:
      For a given (sub)type *T* and gene *g* we calculate the statistical association of mutations in *g* to (sub)type *T* as follows:
        1. Construct the following contingency table:
          .. image:: /_static/subtype_specific_genes_cont_table.*
        2. Perform Fisher's exact test on the above contingency table.
        3. Bonferonni-correct the *p*-value.

    :type gene: string
    :param gene: gene name.
    :type patient2ty: dictionary
    :param patient2ty: mapping of samples to their respective (sub)types (see :func:`load_patient2ty_file`).
    :type ty2numpatients: dictionary
    :param ty2numpatients: mapping of each (sub)type to the number of patients of that subtype.
    :type mutation2patients: dictionary
    :param mutation2patients: mapping of genes to the patients in which they are mutated (see :func:`multi_dendrix.load_mutation_data` for details).

    **Returns:**
      A dictionary that maps each (sub)type to the following tuple:
        (x11, x10, x01, x00, *p*-value, Bonferonni-corrected *p*-value).

    **Examples:**
      A view of the data:
        >>> pnt2ty = {"TCGA-01" : "Luminal A", "TCGA-02" : "Luminal A", "TCGA-03" : "Luminal B"}
        >>> ty2numpatients = {"Luminal A" : 2, "Luminal B" : 1}
        >>> mutation2patients = {"G1" : ["TCGA-01", "TCGA-02", "TCGA-03"],
                                 "G2" : ["TCGA-02"], "G3" : ["TCGA-01", "TCGA-02"]}
      Simple example:
        >>> subtype_specificity("G3", pnt2ty, ty2numpatients, mutation2patients)
        {'Luminal A': [2, 0, 0, 1, 0.33, 0.66], 'Luminal B': [0, 1, 2, 0, 1.0, 1]}

    **See also:** :func:`subtype_analysis`, :func:`ty_contingency_table`.
    """
    tys = sorted(ty2numpatients.keys())
    num_tests = float(len(tys))    

    # Count the number of mutations in the gene in each cancer (sub)type
    ty2mutations = dict([(ty, set()) for ty in tys])
    for sample in mutation2patients[gene]:
        try: ty2mutations[patient2ty[sample]].add( sample )
        except KeyError: continue #  ignore samples with no type

    h = dict()
    for ty in tys:
        cont_table = ty_contingency_table(ty, ty2mutations, tys, ty2numpatients)
        pval = fisher_exact(cont_table)
        corrected_pval = pval * num_tests if pval * num_tests <=1 else 1
        type_mutations, type_normal, non_type_mutations, non_type_normal = cont_table
        
        # Store results
        h[ty] = [ type_mutations, type_normal, non_type_mutations,
                  non_type_normal, pval, corrected_pval ]
    return h

def load_patient2ty_file(patient2ty_file):
    """Loads a file mapping patient IDs to their respective (sub)types.

    :type patient2ty_file: string
    :param patient2ty_file: File location of patient to (sub)type map (see also :doc:`/file_formats`).

    :returns: A tuple mapping patients to their (sub)type, and a list of all (sub)types.

    **Examples:**
      A view of the input data:
        >>> patient_wlst = "patient.wlst"
        >>> open(patient_wlst).readlines()
        ["TCGA-01\\tLuminal A\\n", "TCGA-02\\tLuminal A\\n", "TCGA-03\\tLuminal B\\n"]
      A simple example:
        >>> load_patient2ty_file(patient_wlst)
        {"TCGA-01" : "Luminal A", "TCGA-02" : "Luminal A", "TCGA-03" : "Luminal B"}
    **See also:** :func:`subtype_analysis`, :func:`subtype_specificity`.
    """
    patient2ty = dict([l.rstrip().split("\t")[:2] for l in open(patient2ty_file)
    	              if not l.startswith("#") ])
    tys = sorted(set( [ty for sample, ty in patient2ty.iteritems()]))
    return patient2ty, tys

def keep_significant(gene2specificity, threshold):
    """Removes all associations in the input dictionary that are not significant below a given threshold.

    :type gene2specificity: dictionary
    :param gene2specificity: maps each gene to a dictionary of subtypes, where each subtype is mapped to the result of the statistical test (see output of :func:`subtype_analysis`).
    :type threshold: float
    :param threshold: Significance threshold.

    An example is given in the documentation of :func:`subtype_analysis`.

    **See also:** :func:`subtype_analysis`, :func:`subtype_specificity`.
    """
    sig_gene2specificity = dict()
    for g, ty2analysis in gene2specificity.iteritems():
        h = dict()
        for ty, analysis in ty2analysis.iteritems():
            if analysis[-1] < threshold: # the corrected pvalue
                h[ty] = analysis
		
        if h.keys() != []:
            sig_gene2specificity[g] = h

    return sig_gene2specificity

def create_subtype_tbl(gene2specificity):
    header = ['Gene', '(Sub)Type', 'Type_Mutations', 'Type_Normal',
              'Non_Type_Mutations', 'Non_Type_Normal', 'P_Value',
              'Bonferonni_Corrected_P_Value' ]
    tbl = [  ]
    for gene, ty2analysis in gene2specificity.iteritems():
        for ty, analysis in ty2analysis.iteritems():
            tbl.append( [gene, ty] + analysis )
    
    # Sort rows by (sub)type, then p-value, then gene name
    tbl.sort(key=lambda arr: (arr[1], arr[-1], arr[0]))
    tbl = [header] + tbl
    return [ map(str, row) for row in tbl ]

def subtype_analysis(mutation_data, patient_whitelist, threshold=1.0):
    """Performs analysis for subtype-specific genes or mutation classes in given mutation data.

    See :func:`subtype_specificity` for details of the statistical test.

    :type mutation_data: tuple
    :param mutation_data: Mutation data as output from :func:`multi_dendrix.load_mutation_data`.
    :type patient_whitelist: string
    :param patient_whitelist: Location of patient whitelist file that also includes the (sub)type of each patient (see :doc:`/file_formats` for details).
    :type threshold: float
    :param threshold: Significance threshold for *Bonferonni-corrected* *p*-values (default: 1.0).

    **Returns:**
      Mapping of genes (mutation classes) to the following tuple for each (sub)type  *T*:
        * The number of patients of type *T* with mutations in the gene. 
        * The number of patients of type *T* *without* mutations in the gene.
        * The number of patients *not* of type *T* with mutations in the gene.
        * The number of patietns *not* of type *T* *without* mutations in the gene.
        * The (uncorrected) Fisher's exact test *p*-value of the association of the gene with *T*.
        * The Bonferonni-corrected *p*-value.

    **Examples:**
      A view of the data:
        >>> mutation_data = (2, 3, ["G1", "G2", "G3"], ["TCGA-01", "TCGA-02", "TCGA-03"],
            {"G1" : ["TCGA-01", "TCGA-02", "TCGA-03"], "G2" : ["TCGA-02"],
             "G3" : ["TCGA-01", "TCGA-02"]},
            {"TCGA-01" : ["G1", "G3"], "TCGA-02" : ["G1", "G2", "G3"], "TCGA-03" : ["G1"]})
        >>> patient_wlst = "patient.wlst"
        >>> open(patient_wlst).readlines()
        ["TCGA-01\\tLuminal A\\n", "TCGA-02\\tLuminal A\\n", "TCGA-03\\tLuminal B\\n"]
      Example with no significance threshold:
        >>> subtype_analysis(mutation_data, patient_wlst)
        {'G3': {'Luminal A': [2, 0, 0, 1, 0.33, 0.66], 'Luminal B': [0, 1, 2, 0, 1.0, 1]},
        'G2': {'Luminal A': [1, 1, 0, 1, 0.66, 1], 'Luminal B': [0, 1, 1, 1, 1.0, 1]},
        'G1': {'Luminal A': [2, 0, 1, 0, 1.0, 1], 'Luminal B': [1, 0, 2, 0, 1.0, 1]}}
      Example with a (unusually high) significance threshold of 0.7:
        >>> subtype_analysis(mutation_data, patient_wlst, 0.7)
        {'G3': {'Luminal A': [2, 0, 0, 1, 0.33, 0.66]}}

    **See also:** :func:`subtype_specificity`, :func:`keep_significant`.

    """
	# Parse mutation data and load patient2ty file
    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data
    patient2ty, tys = load_patient2ty_file(patient_whitelist)

    # Count the number of samples from each cancer (sub)type
    ty2numpatients = dict([(ty, 0) for ty in tys])
    for sample, ty in patient2ty.iteritems(): ty2numpatients[ty] += 1

    gene_specificity = [ subtype_specificity(g, patient2ty, ty2numpatients,
                                              mutation2patients)
                         for g in genes ]
    gene2specificity = dict(zip(genes, gene_specificity))

    # Prune list if required
    if threshold   < 1.0:
        gene2specificity = keep_significant(gene2specificity, threshold)

    return gene2specificity


def run(args):
    """Analyze the given mutation data for subtype-specific mutations, and output the results to file."""
	# Load mutation data
    from .. import multi_dendrix as Multi
    
    include = Multi.white_and_blacklisting(args.patient_whitelist,
     args.patient_blacklist, args.gene_whitelist, args.gene_blacklist)
    gene2include, patient2include = include

    mutation_data = Multi.load_mutation_data_w_cutoff(args.mutation_matrix,
        patient2include, gene2include, args.cutoff)

    # Conduct subtype analysis
    threshold = None if args.all else args.sig_threshold
    gene2specificity = subtype_analysis(mutation_data, args.patient_whitelist,
    	                                threshold)

    # Create TSV table to output results
    subtype_tbl = create_subtype_tbl(gene2specificity)
    subtype_output = "\n".join([ "\t".join(row) for row in subtype_tbl ])

    # Output results to file
    if args.output_file:
        open(args.output_file, 'w').write( subtype_output )
    else:
        print subtype_output

if __name__ == "__main__": run(parse_args())