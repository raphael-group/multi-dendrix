#!/bin/sh

# ##############################################################################
# Script to run the Multi-Dendrix pipeline on the example GBM(2008) data.
# To get started, run this script to make sure that it works on your machine.
# Then, I would copy this script and use it as a template for your own scripts.
#
# For more examples and information, check out the documentation!
# http://compbio.cs.brown.edu/projects/multi-dendrix/docs
################################################################################


# General parameters
OUTPUT_DIR=output/GBM
VERBOSE="-v"

# Parameters for the  Multi-Dendrix ILP
KMIN=3
KMAX=5
TMIN=2
TMAX=4
NAME=GBM
MUTATION_MATRIX=examples/mutation_data/GBM_2008/GBM_2008.m2
GENE_BLACKLIST=examples/mutation_data/fishy_genes.glst

# Parameters related to the evaluation
PERMUTE="--permute"
PPI_NETWORK=examples/network/iref_edgelist
NUM_PERMUTED_MATRICES=1

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

./multi_dendrix_pipeline.py -o $OUTPUT_DIR $VERBOSE -k_min $KMIN -k_max $KMAX \
                            -t_min $TMIN -t_max $TMAX -n $NAME \
                            -m $MUTATION_MATRIX -bg $GENE_BLACKLIST $PERMUTE\
                            -ppi $PPI_NETWORK \
                            --num_permuted_matrices $NUM_PERMUTED_MATRICES

