Multi-Dendrix
=============

This Python package implements the Multi-Dendrix algorithm for analyzing driver
pathways in cancer mutation data ([*PLoS Comp Bio, 2013*](http://goo.gl/fF0N4)).
It was developed by the [Raphael research group](http://compbio.cs.brown.edu)
in the Center for Computational Molecular Biology at Brown University.

*Update: 2014/06/28*

This package also includes a Python package that implements a Markov chain Monte
Carlo (MCMC) version of Multi-Dendrix. Multi-Dendrix MCMC does not require IBM's
CPLEX, and is thus completely open-source to use. We include a brief description
of how to run Multi-Dendrix MCMC at the end of this README.

Install
--------
After pulling from Github, simply run:

  `python setup.py install`

System requirements
--------------------
Multi-Dendrix was developed on 64-bit Debian Linux. It has not been tested
on other systems.

The Multi-Dendrix package requires the following Python modules:
  * [IBM's CPLEX] (http://goo.gl/dJV6f).
  * [NetworkX] (http://networkx.github.com).
  * Either SciPy >= version 0.11 or [fisher0.1.4](http://goo.gl/zYrLr).

In addition, for web output with network figures, Multi-Dendrix requires the
installation of [GraphViz](http://www.graphviz.org/).

Usage
------
To run the full Multi-Dendrix pipeline on the provided data,
   
   `./run_multi_dendrix.sh`

For other uses, please refer to our documentation (see below).

Documentation
--------------
We offer extensive documentation for running Multi-Dendrix, as well as creating
custom analysis pipelines, at
[http://mdml.github.com/multi-dendrix](http://mdml.github.com/multi-dendrix).

Support
-----------------
Please visit our [Google Group](https://groups.google.com/forum/#!forum/dendrix)
to post questions and view discussions from other users.

Multi-Dendrix MCMC
----------------------
The Multi-Dendrix MCMC algorithm samples collections of *t* gene sets of size *k* (
note that *k* is fixed, unlike in the Multi-Dendrix ILP) in proportion to their weight. 
The Multi-Dendrix release now includes a Python package where the MCMC algorithm is implemented as a Python C extension. To run Multi-Dendrix MCMC:

1. First, compile the C code to build the Multi-Dendrix MCMC Python module:

       cd multi_dendrix_mcmc
       python compile.py build
2. Then, get a complete list of arguments to the Multi-Dendrix MCMC program by
   running:
   
       cd ../
       python run_multi_dendrix_mcmc.py -h
       
The output of Multi-Dendrix MCMC is a tab-separated file that lists the collections
sampled by the MCMC algorithm. The first line of each collection lists the sampling
frequency, weight *W'*, permutation P-value (if computed), and the first gene set
in the collection and its weight. The following *t-1* lines list the remaining gene
sets and their weights.

Multi-Dendrix MCMC computes the permutation *P*-values for the *i*th highest weight
collection using the same procedure for permuting mutation data described in the
Multi-Dendrix paper. The MCMC algorithm is run on each of these permuted datasets,
and the weight of the *i*th highest weight collection from real data is compared
to the weight of the *i*th highest weight collection across the permtued datasets.
**Multi-Dendrix MCMC only computes the permutation P-value for the top 25 highest
weight collections from a given run.**

We thank [Troy D. Hanson](http://troydhanson.github.io/) for the
[uthash library](http://troydhanson.github.io/uthash/) used by the Multi-Dendrix
MCMC C extension.