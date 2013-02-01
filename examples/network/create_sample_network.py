# Load required modules
MULTI_DIR="/home/mdml/Desktop/research/dendrix/multi-dendrix"
import sys, random
sys.path.insert(1, MULTI_DIR)
from evaluate import *
import networkx as nx

# Load iRef
iRef = load_PPI("iRef")

# Choose genes to keep
NUM_GENES=500
genes = random.sample(iRef.nodes(), NUM_GENES)

# Induce subnetwork and output edgelist to file
G = nx.subgraph(iRef, genes)
nx.write_edgelist(iRef, 'iref_edgelist')
