Multi-Dendrix
=============

This Python package implements the Multi-Dendrix algorithm for analyzing driver pathways in cancer mutation data (paper in submission). It was developed by the Raphael research group in the Center for Computational Molecular Biology at Brown University.

Install
--------
After pulling from Github, simply run:
  >>> python setup.py install

System requirements
--------------------
Multi-Dendrix was developed on 64-bit Debian Linux. It has not been tested on other systems.

The Multi-Dendrix package requires the following Python modules:
  * `IBM's CPLEX <http://goo.gl/dJV6f>`_.
  * `NetworkX <networkx.github.com>`_.
  * Either SciPy v0.11 or `fisher0.1.4 <http://goo.gl/zYrLr>`_.

In addition, for web output with network figures, Multi-Dendrix requires the installation of `GraphViz <http://www.graphviz.org/>`_.

Usage
------
To run the full Multi-Dendrix pipeline on the provided data,
   >>> ./run_multi_dendrix.sh
For other uses, please refer to our documentation (see below).

Documentation
--------------
We offer extensive documentation for running Multi-Dendrix, as well as creating custom analysis pipelines, at http://multi-dendrix.github.com.

