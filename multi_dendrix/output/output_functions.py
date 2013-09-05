#!/usr/bin/python

# Load required modules
import textwrap, os
from time import strftime
from networkx import connected_components, write_dot

###############################################################################
# Functions for formatting results in Python tables

def format_collection_tbl(args, collection, weights, set_eval):
    """Creates a table from a collection of gene sets of the following form:
    ======   ============   ===================  ====================
    Weight   Network Test   Permutation P-Value  Genes
    (*int*)  (*float*)      (*float*)            (*list of strings*)
    ======   ============   ===================  ====================
    """
    if set_eval: test, stat, pval, gene_set_results = set_eval
    else: test = 'Permutation_Test'
    
    tbl = [ ['Weight', test, 'Permutation_P_Value', 'Gene Set'] ]
    for i in range(len(collection)):
        if set_eval: _, stat, pval = gene_set_results[i]
        else: stat, pval = "--", "--"
        tbl.append( [ weights[i], stat, pval, ', '.join(collection[i])] )
    
    return tbl

def create_collection_tbls(args, collections, core_modules, evaluation):
    """Creates a table using :func:`format_collection_tbl` for
    each collection."""
    # Check if an evaluation was performed
    network_results, matrix_results = evaluation

    tables = dict( [ (t, {}) for t in collections.keys() ] )
    for t in collections.keys():
        for k_max in collections[t].keys():
            # For each gene set, output the weight, the gene set,
            # and its evaluation
            collection, weights = collections[t][k_max]
            if args.network_test: set_eval = network_results[t][k_max]
            else: set_eval = None
            tables[t][k_max] = format_collection_tbl(args, collection, weights,
                                                      set_eval)

    # Add the core module information
    if args.network_test: set_eval = network_results["core_modules"]
    else: set_eval = None
            
    weights = [ "--" for i in range(len(core_modules)) ]
    tables["core_modules"] = format_collection_tbl(args, core_modules,
                                                   weights, set_eval)

    return tables

def create_matrix_results_tbl(matrix_results):
    """Creates a table from the results of the matrix permutation tests
    of the following form:
    ======   =======  ===================
    *t*      *kmax*   Permutation P-Value
    (*int*)  (*int*)  (*float*)          
    ======   =======  ===================
    """
    tbl = [ ['t', 'k_max', 'Permutation_P_value']]
    for t in matrix_results.keys():
        for k_max in matrix_results[t].keys():
            tbl.append( [ t, k_max, matrix_results[t][k_max] ] )
    return [map(str, row) for row in tbl]

def create_network_results_tbl(network_results):
    """Creates a table from the results of the network test (either *Direct 
    Interactions* or *Average Pairwise Distance* test of the following form:
    ======   =======  =======================  =====================
    *t*      *kmax*   Network Test Statistic   Network Test P-Value
    (*int*)  (*int*)  (*float*)                (*float*)
    ======   =======  =======================  =====================
    """
    tbl = [ ['t', 'k_max', 'Test', 'P_value'] ]
    for t in network_results.keys():
        if t == "core_modules":
            test_name, statistic, pval, gene_set_results = network_results[t]
            tbl.append( ['Core_Modules', '', statistic, pval] )
            continue

        for k_max in network_results[t].keys():
           results = network_results[t][k_max]
           test_name, statistic, pval, gene_set_results = results
           tbl.append( [ t, k_max, statistic, pval ] )

        if tbl[0][2] == 'Test': tbl[0][2] = test_name

    return [map(str, row) for row in tbl]

def create_params_tbl(args, mutation_data):
    """Creates a table from the results of the matrix permutation tests
    of the following form:
    ==========  ===================================
    Parameter   Value
    (*str*)      (*str* / *float* / *int* / *list*)
    ==========  ===================================
    """
    # Create arrays of parameters that will be required later
    ts = range(args.min_num_gene_sets, args.max_num_gene_sets+1)
    ks = range(args.min_gene_set_size, args.max_gene_set_size+1)

    m, n, genes, patients, mutation2patients, patient2mutations = mutation_data

    tbl = [ ['Param', 'Value'] ]
    tbl.append( [ 'General parameters' ] )
    tbl.append( [ 'Output directory', args.output_dir ])
    tbl.append( [ 'Input data' ])
    tbl.append( [ 'Dataset', args.db_name ])
    tbl.append( [ 'Mutation matrix', args.mutation_matrix ])
    tbl.append( [ 'Minimum mutation frequency', args.cutoff ])
    tbl.append( [ 'Number of patients', n ])
    tbl.append( [ 'Number of genes', m ])
    tbl.append( [ 'Patient whitelist', args.patient_whitelist ])
    tbl.append( [ 'Patient blacklist', args.patient_blacklist ])
    tbl.append( [ 'Gene whitelist', args.gene_whitelist ])
    tbl.append( [ 'Gene blacklist', args.gene_blacklist ])

    tbl.append( [ 'ILP Parameters' ] )
    tbl.append( [ 'Range of number of gene sets', ', '.join(map(str, ts)) ])
    tbl.append( [ 'Range of max. gene set sizes', ', '.join(map(str, ks)) ])
    tbl.append( [ 'k<sub>min</sub>', args.min_num_gene_sets ])
    tbl.append( [ 'k<sub>max</sub>', args.max_num_gene_sets ])
    tbl.append( [ 'Exclusivity weight (&alpha;)', args.alpha ])
    tbl.append( [ 'Number of overlaps between gene sets (&delta;)',
                  args.delta ] )
    tbl.append( [ 'Number of gene sets per gene (&lambda;)',
                  args.lmbda ])
    tbl.append( [ 'Stability threshold of core modules',
                   args.stability_threshold ] )

    if args.subtypes:
        tbl.append( [ 'Subtype analysis' ] )
        tbl.append( [ 'Subtype significance threshold',
                      args.subtype_sig_threshold ])

    if args.network_test or args.weight_test:
        tbl.append( [ 'Evaluation parameters' ] )

    if args.network_test:
        tbl.append( [ 'Network location', args.network_edgelist ] )
        tbl.append( [ 'Number of permuted networks',
                      args.num_permuted_networks ] )
        tbl.append( [ 'Permuted networks directory',
                      args.permuted_networks_dir ] )
        tbl.append( [ 'Edge swap parameter Q', args.Q ] )
        if args.distance:
            tbl.append( [ 'Network test statistic', 'Avg. pairwise distance' ])
        else:
            tbl.append( [ 'Network test statistic', 'No. of interactions'])

    if args.weight_test:
        tbl.append( [ 'Number of permuted matrices',
                      args.num_permuted_matrices ] )
        tbl.append( [ 'Permuted matrices directory',
                      args.permuted_matrices_dir ] )

    return [ map(str, row) for row in tbl ]

###############################################################################
# Helper functions for outputting web pages

def page_header(title, today=strftime("%a, %d %b %Y %H:%M:%S")):
    """Return a custom header for an HTML page. Uses Twitter Bootstrap."""
    return textwrap.dedent(
        '''\
        <html>
        <head>
        <meta http-equiv="content-type" content="text/html;charset=UTF-8" />
        <link href="http://netdna.bootstrapcdn.com/twitter-bootstrap/2.2.1/css
        /bootstrap-combined.min.css" rel="stylesheet">
        <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/jquery
        .min.js"></script>
        <script src="http://netdna.bootstrapcdn.com/twitter-bootstrap/2.2.1/
        js/bootstrap.min.js"></script>
        <style>
        .body{ padding:5px; }
        #date{ float:right; padding-top:10px; }
        .container{ padding:5px; }
        #title{ margin-left:5px; }

        a.back-link{ float:left; }
        div.stats-links{ float:right; }

        table#params td{ text-align:left; }
        table#subtypes td, table#subtypes th{ text-align:center; }

        th.params-subhead{ text-align:left; text-decoration:underline; }
        </style>
        <title>%(title)s.</title>
        </head>
        <body>
        <div class="navbar navbar-inverse navbar-fixed-top" id="navigation">
        <div class="navbar-inner">
        <span id="date">%(today)s</span>
        <div>
        <div class="brand" id="title">%(title)s.</div>
        </div>
        </div>
        </div><br/><br/><div class="body">
        ''' % { "today" : today, "title" : title })

def format_pval(pval, num_random_samples):
    """Formats the *p*-value of a permutation test.
    If *pval* is not zero, returns it with a prefix.
    Otherwise, it returns *P* < 1/ *num_random_samples*."""
    if pval != 0: return "<i>P</i> = " + format(pval, 'g')
    else: return "<i>P</i> < " + format(1. / num_random_samples, 'g')

###############################################################################
# Functions for outputting the individual pages in the Multi-Dendrix website

def params_page(params_tbl):
    """Function for piecing together the web page that lists the parameters
    of the run."""
    page = page_header('Multi-Dendrix Params')
    page += "<a href='./index.html' class='back-link'>Back</a>\n<br/>"
    page += "<h1>Parameters</h1>\n"
    page += "<table class='table table-striped table-condensed' id='params'>\n"
    page += "<tr>" + "".join( ["<th>" + th + "</th>" for th in params_tbl[0]])
    for row in params_tbl[1:]:
        if len(row) == 1:
            page += '<tr><th class="params-subhead" colspan="2">' + row[0]
        else:
            page += "<tr><th>" + row[0] + "</th><td>" + row[1]
        page += "</td></tr>\n"
    page += "</table>\n"
    return page + "</body></html>"

def subtypes_page(args, subtype_tbl):
    """Function for piecing together the web page that lists subtype-specifc
    genes / mutation classes."""
    page = page_header('Subtype-specific genes')
    page += "<div class='stats-links'>\n<a href='./core_modules."\
            "html'>Core modules</a></div>"
    page += "<a href='./index.html' class='back-link'>Back</a>\n<br/>"
    page += "<h1>Subtype-specific genes</h1>\n"
    page += "Significance threshold: <i>p</i> < "
    page += str(args.subtype_sig_threshold) + "\n"
    page += "<table class='table table-striped table-condensed' id='subtypes'>"
    page += "<tr>" + "".join( ["<th>" + th + "</th>" for th in subtype_tbl[0]])
    for row in subtype_tbl[1:]:
        page += "<tr><th>" + row[0] + "</th>"
        page += "".join( [ "<td>" + cell + "</td>" for cell in row[1:] ] )
        page += "</tr>\n"
    page += "</table>\n"
    return page + "</body></html>"

def collection_page(args, t, k_max, collection, weights, net_eval, mtx_pval):
    """Function for piecing together the web page that dislpays a single
    collection (and its evaluation)."""
    page = page_header(args.db_name + " t=" + str(t) + ", kmax=" + str(k_max))
    page += "<a href='../index.html' class='back-link'>Back</a>\n<br/>"
    
    # Output the collection and their evaluation
    page += "<h1>Collection</h1>\n"
    page += "<table class='table table-striped table-condensed'>"
    page += "<tr><th>Genes</th><th>Weight</th>"
    if net_eval:
        test, collect_stat, collect_pval, gene_set_results = net_eval
        page += "<th>" + test + "</th><th><i>P</i>-value</th>"
    page += "</tr>\n"
    for i in range(len(collection)):
        page += "<tr><td>" + ', '.join(collection[i]) + "</td>"
        page += "<td>" + str(weights[i]) + "</td>"
        if args.network_test:
            _, set_stat, set_pval = gene_set_results[i]
            page += "<td>" + format(set_stat, 'g') + "</td>"
            page += "<td>" + format_pval(float(set_pval), args.num_permuted_networks)
            page += "</td>"
        page += "</tr>\n"
    page += "\n</table>"

    # Evaluation (if necessary)
    if net_eval:
        page += "<h1>Evaluation of collection</h1>\n"
        page += "<table class='table table-striped table-condensed'>\n"
        page += "<tr><th>Test</th><th>Statistic</th><th><i>P</i>-value</th>"\
                "</tr>\n"
        page += "<tr><td><i>Network permutation test</i></td><td>"
        page += format(collect_stat, 'g') + "</td><td>"
        page += format_pval(collect_pval, args.num_permuted_networks)
        page += "</td><tr>\n"
        
    if args.weight_test:
        page += "<tr><td><i>Matrix permutation test</i></td>"
        page += "<td>--</td><td>" + format_pval(mtx_pval, args.num_permuted_matrices)
        page += "</td></tr>\n"
    
    return page + "</table>\n</body></html>"

def index_page(args, collections, runtime, (network_eval, matrix_eval)):
    """Function for piecing together the web page that lists all the collections
    produced by a run, as well as links to all the individual pages."""
    page = page_header("Multi-Dendrix Results on " + args.db_name)
    page += "<a href='./core_modules.html' class='back-link'>Core "\
            "Modules</a>\n"
    page += "<div class='stats-links'>\n"
    if args.subtypes:
        page += " <a href='./subtype_analysis.html'>Subtype analysis</a> | "
    page += "<a href='./params.html' class='stats-link'>Params</a>\n</div>\n"

    # Output collections
    page += "<br/>\n<h1>Collections</h1>\n"
    page += "<b>Runtime: </b>" + '%.3f' % runtime + " seconds."
    if not (args.weight_test or args.network_test):
        page += "<br/><i>Evaluation not performed.</i>"
    page  += "<table class='table table-striped table-condensed'>\n"
    page += "<tr><th><i>t</i></th><th><i>k</i><sub>max</sub></th><th>"\
            "Gene Sets</th><th>Weights</th>"
    if args.network_test:
        page += "<th>Network permutation test</th>"
    if args.weight_test:
        page += "<th>Matrix permutation test</th>"
    page += "<th>Link</th></tr>\n"
    for t in [ key for key in collections.keys() if key != "core_modules" ]:
        for k_max in collections[t].keys():
            if args.network_test:
                _, set_stat, set_pval, _ = network_eval[t][k_max]
                net_pval = format_pval(set_pval, args.num_permuted_networks)
            if args.weight_test:
                mtx_pval = format_pval(matrix_eval[t][k_max],
                                       args.num_permuted_matrices)

            collection, weights = collections[t][k_max]

            page += "<tr><td>" + str(t) + "</td><td>" + str(k_max) + "</td>"
            page += "<td>" + "<br/>".join( [", ".join(p) for p in collection] )
            page += "</td><td>" + "<br/>".join(map(str, weights)) + "</td>"
            if args.network_test:
                page += "<td>" + net_pval + "</td>"
            if args.weight_test:
                page += "<td>" + mtx_pval + "</td>"

            run_name = "-".join([args.db_name, "t"+str(t), "kmax"+str(k_max)])
            run_link = "<a href='collections/" + run_name + ".html'>>></a>"
            page += "<td>" + run_link + "</td></tr>\n"

    return page + "</body></html>"

###############################################################################
# Functions for drawing the core module diagram and creating the core module 
# page

def style_nodes(nodes, subtype_genes):
    """Changes the style of NetworkX Nodes based on whether genes are enriched for subtype-specific mutations (or not).
    Genes that are enriched for subtype-specific mutations are rectangles.
    Genes that are *not* enriched for subtype-specific mutations are rounded rectangles.
    """
    # Change genes' shapes based on whether their mutations are enriched 
    # for a particular subtype
    subtype_specific = dict(style="filled")
    general = dict(style="rounded, filled")

    general_nodes = [ (n, general) for n in nodes if n not in subtype_genes ]
    subtype_nodes = [ (n, subtype_specific) for n in subtype_genes ]

    return general_nodes + subtype_nodes

def style_edges(edges, subtype_genes, gene2specificity):
    """Changes the style of NetworkX Edges based on whether a pair of genes are enriched for *different* subtype-specific mutations.
    Edges connecting genes that are not enriched for subtype-specific mutations, OR are enriched for the same subtype, are connected by solid lines.
    Edges connecting genes that are enriched for subtype-specific mutations of different subtypes are dashed.
    """
    if gene2specificity:
        gene2tys = dict([(g, set(gene2specificity[g].keys()))
                         for g in subtype_genes])

    styled_edges = [ ]
    for u, v, h in edges:
        style = dict(penwidth=h["weight"], color="blue", label=h["weight"])
        if u in subtype_genes and v in subtype_genes:
            if len(gene2tys[u] | gene2tys[v]) > 1:
                style["style"] = "dashed"
        styled_edges.append( (u, v, style) )

    return styled_edges

def core_module_drawing(G, gene2specificity, dot_fh="core_modules", ext="svg"):
    """Draws a network representation of the core modules using GraphViz.
    For information on how the core modules are identified, see :func:`core_modules.extract_core_modules`.
    """
    H = G.copy()
    H.add_node("graph", bgcolor="transparent")
    
    if gene2specificity:
        subtype_genes = set( gene2specificity.keys() ) & set( G.nodes() )
    else: subtype_genes = []

    H.add_nodes_from(style_nodes(G.nodes(), subtype_genes), shape="rect",
                     fontcolor="black", fontsize=11)
    H.add_edges_from(style_edges(G.edges(data=True), subtype_genes,
                                 gene2specificity))

    H.remove_nodes_from([n for n in H.nodes() if len(H[n]) == 0 ])

    write_dot(H, dot_fh + ".dot")
    os.system("dot -T" + ext + " " + dot_fh + ".dot -o " + dot_fh  + "." + ext)
    svg = open(dot_fh + "." + ext).read()
    os.remove(dot_fh + "." + ext)

    return svg

def core_module_page(args, G, output_dir, net_eval, gene2specificity):
    """Function for piecing together the various elements that make up the
    web page summarizing the core modules found for a given run."""
    page = page_header(args.db_name + " Multi-Dendrix Core Modules")
    if gene2specificity:
        page += "<div class='stats-links'>\n<a href='./subtype_analysis."\
                "html'>Subtype analysis</a></div>"
    page += "<a href='./index.html' class='back-link'>Back</a>\n<br/>"
    # 
    components = connected_components(G)
    modules = [ subG for subG in components if len(subG) > 2 ]
    singletons = [ subG[0] for subG in components if len(subG) == 1 ]

    # Output the modules as a graph
    page += "<h1>Module visualization</h1>\n"
    page += core_module_drawing(G, gene2specificity,
                                output_dir + "/core_modules")

    # Add the modules as rows in a table below
    page += "<hr/><h1>Core modules</h1>\n"
    page += "<table class='table table-striped table-condensed'>"
    page += "<tr><th>Genes</th>"
    if net_eval:
        test, set_stat, set_pval, module_results = net_eval
        page += "<th>" + test + "</th><th><i>P</i>-value</th>"
    page += "</tr>\n"

    for i in range(len(modules)):
        page += "<tr><td>" + ', '.join(modules[i]) + "</td>"
        if net_eval:
            _, path_stat, path_pval = module_results[i]
            page += "<td>" + format(path_stat, 'g') + "</td><td>"
            page += format_pval(path_pval, args.num_permuted_networks)
            page += "</td>"
        page += "</tr>\n"
    page += "\n</table>"
    if len(singletons) != 0:
        page += "<b>Modules of size 1:</b> " + ", ".join(singletons)

    # Evaluation (if necessary)
    if net_eval:
        page += "<h1>Evaluation of module set</h1>\n"
        page += "<table class='table table-striped table-condensed'>\n"
        page += "<tr><th>Test</th><th>Statistic</th><th><i>P</i>-value</th>"\
                "</tr>\n<tr><td><i>"
        if args.distance: page += "Average Pairwise Distance"
        else: page += "Direct Interactions"
        page += " Test</i></td><td>" + format(set_stat, 'g') + "</td><td>"
        page += format_pval(set_pval, args.num_permuted_networks)
        page += "</td></tr>\n</table>\n"
    

    return page + "</body></html>"

###############################################################################
# Main functions for outputting results to text or html
def output_to_text(args, collection_tbls, runtime, params_tbl, network_tbl,
                   matrix_tbl, subtype_tbl):
    """Outputs the results of a Multi-Dendrix run to text files."""
    # Make sure output directory exists
    output_dir = args.output_dir + "/txt/"
    os.system('mkdir -p ' + output_dir + "/collections")

    # For each parameter set, output the gene sets and their evaluations
    for t in collection_tbls.keys():
        if t == "core_modules":
            rows = [ "\t".join(map(str, row)) for row in collection_tbls[t]]
            output_file = output_dir + "/collections/core_modules.tsv"
            open(output_file, 'w').write( "\n".join( rows ))
            continue

        for k_max in collection_tbls[t].keys():
            rows = [ "\t".join(map(str, row))
                     for row in collection_tbls[t][k_max] ]

            run_name = "-".join([args.db_name, "t" + str(t), "kmax" + str(k_max)])
            output_file = output_dir + "collections/" + run_name + '.tsv'

            open(output_file, 'w').write( '\n'.join(rows) )

    # Output params file
    params_output = "\n".join( [ "\t".join(row) for row in params_tbl ] )
    open(output_dir + "params.tsv", "w").write( params_output )

    # Output Multi-Dendrix runtime as single file
    open(output_dir + "runtime", "w").write( str(runtime) )

    # Output network permutation test results
    if network_tbl:
        network_output = "\n".join( [ "\t".join(row) for row in network_tbl ])
        open(output_dir + "network_permutation_test.tsv", "w").write(
            network_output)

    # Output matrix permutation test results
    if matrix_tbl:
        matrix_output = "\n".join( [ "\t".join(row) for row in matrix_tbl ])
        open(output_dir + "matrix_permutation_test.tsv", "w").write(
            matrix_output)

    # Output subtype analysis
    if subtype_tbl:
        subtype_output = "\n".join( [ "\t".join(row) for row in subtype_tbl ])
        open(output_dir + "subtype_specific_genes.tsv", "w").write(
            subtype_output)

def output_to_html(args, collections, runtime, core_module_graph, evaluation,
                   params_tbl, subtype_tbl, gene2specificity):
    """Outputs the results of a Multi-Dendrix run as a web site."""
    network_eval, matrix_eval = evaluation

    # Make sure output directory exists
    output_dir = args.output_dir + "/html/"
    os.system("mkdir -p " + output_dir + "collections")

    # Output index page
    open(output_dir + "/index.html", "w").write(
        index_page(args, collections, runtime, evaluation)
        )

    # Output core modules
    if args.network_test: core_module_eval = network_eval["core_modules"]
    else: core_module_eval = None
    open(output_dir + "/core_modules.html", "w").write(
        core_module_page(args, core_module_graph, output_dir,
                         core_module_eval, gene2specificity))

    # Output each collection
    for t in [ key for key in collections.keys() if key != "core_modules" ]:
        for k_max in collections[t].keys():
            if args.network_test: net_eval = network_eval[t][k_max]
            else: net_eval = None

            if args.weight_test: mtx_pval = matrix_eval[t][k_max]
            else: mtx_pval = None

            collection, weights = collections[t][k_max]

            run_name = "-".join([args.db_name, "t" + str(t), "kmax" + str(k_max)])
            open(output_dir + "collections/" + run_name + ".html", "w").write(
                collection_page(args, t, k_max, collection, weights,
                                 net_eval, mtx_pval)
                )

    # Output params page
    open(output_dir + "params.html", "w").write( params_page(params_tbl) )

    # Output subtypes table (if necessary)
    if subtype_tbl:
        open(output_dir + "subtype_analysis.html", "w").write(
            subtypes_page(args, subtype_tbl) )
