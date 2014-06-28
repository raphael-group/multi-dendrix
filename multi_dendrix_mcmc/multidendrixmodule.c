#include <Python.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include "uthash.h"

/*****************************************************************************/
/* Utility array functions */
// Returns the index of n in array where k is the length of array, or -1 if n is
// not in array.
int indexOf(int n, int *array, int k){
  int i;
  for (i=0; i < k; i++) if(array[i] == n) return i;
  return -1;
}

// Copy an array of doubles
void copyArr(double *A, double *B, int n){
  int i;
  for (i = 0; i < n; i++) B[i] = A[i];
}

// Sum an array of doubles
double sum(double *arr, int n){
  int i;
  double total = 0;
  for (i = 0; i < n; i++) total += arr[i];
  return total;
}

/*****************************************************************************/
/* Computing the weight */
// Structs for the hash
typedef struct{
    int genes[10];
} geneset_t;

typedef struct {
    geneset_t id; /* REQUIRED: the key */
    double weight;
    UT_hash_handle hh; /* REQUIRED: makes this structure hashable */
} weight_t;

// Global store of the weights of seen gene sets
weight_t *weightHash = NULL;
int lookups = 0;
int calculations = 0;

// Dendrix weight function
double W(int start, int k, int n, bool **A, int *set){
  // Declarations
  double set_coverage, total_coverage, sample_mutations;
  int i, j;

  // Iterate through the samples, 
  set_coverage   = 0.0;
  total_coverage = 0.0; 

  for(j = 0; j < n; j++){
    sample_mutations = 0.0;
    for(i = start; i < start + k; i++){
      if (A[set[i]][j]) sample_mutations += 1.0;
    }
    total_coverage += sample_mutations;
    if (sample_mutations != 0.0) set_coverage += 1.0;
  }
  
  return 2.0 * set_coverage - total_coverage;
}

// Sort a pair of integers in ascending order
int ascending(const void * elem1, const void * elem2){
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

// Get the weight for a geneset, either from the hash
// or by computing it then adding it to the hash
double getWeight(int *set, int k, int start, int n, bool **A){
  // Declarations
  geneset_t *key;
  weight_t *w;
  int i;

  // Generate the key by sorting the arr
  key = malloc( sizeof(geneset_t) );
  memset(key, 0, sizeof(geneset_t));
  for (i = 0; i < k; i++) key->genes[i] = set[i + start];
  qsort(key->genes, k, sizeof(int), ascending);

  // Try and find the geneset
  HASH_FIND(hh, weightHash, key, sizeof(geneset_t), w);  /* id already in the hash? */

  // If it's not there, add it
  if (w == NULL) {
    calculations += 1;
    w = malloc(sizeof(weight_t));
    w->id = *key;
    w->weight = W(start, k, n, A, set);
    HASH_ADD( hh, weightHash, id, sizeof(geneset_t), w );  /* id: name of key field */
  }
  else{
    lookups += 1;
  }

  // Return the gene set's weight
  return w->weight;
}

/*****************************************************************************/
// Generate a collection of t random sets of k genes (integers)
void rand_soln(int *set, double *weights, int m, int n, int t, int k, bool **A){
  // Declarations
  int i, r;

  // Select t random gene sets of size k, making sure
  // there are no duplicates
  for(i = 0; i < t * k; i++){
    do{
      r = rand() % m;
    } while(indexOf(r, set, i) != -1);
    set[i] = r;
    if ((i+1) % k == 0){
      weights[i/k] = getWeight(set, k, i/k*k, n, A);  
    }
  }
}

// Multi-Dendrix MCMC
void mcmc(bool **A, int m, int n, int k, int t, int num_iters, int step_len,
		  int **gene_sets, double **set_weights, int verbose){
  // Declarations
  int i, j, step_index, next_gene, gene_index, num_genes = k * t;
  int other_set_to_exchange, set_to_exchange, to_exchange, exchanged;
  int steps, prog_step = num_iters/72; // for progress bar
  int num_swapped = 0, num_within_swap = 0; // MCMC stats
  double log_acceptance_prob, c, logcoin;
  double *weights, *next_weights, next_weight, current_weight = 0.0;
  int *set;
  clock_t begin = clock(), end;

  // Allocations
  set          = malloc(sizeof(int) * num_genes);
  weights      = malloc(sizeof(weights[0]) * t);
  next_weights = malloc(sizeof(next_weights[0]) * t);

  // Intialize random number generator with the current time
  time_t T;
  srand((unsigned) time(&T));
  
  // Initialize the Markov chain
  rand_soln(set, weights, m, n, t, k, A);
  current_weight = sum(weights, t);

  // Run MCMC process
  c = 0.5;
  step_index = 0;
  for (i = 0; i < num_iters; i++){
    // Determine the genes to swap in / out of the solution
    next_gene        = (double) (rand() % m);
    to_exchange      = rand() % num_genes;
    set_to_exchange  = to_exchange / k;
    gene_index       = indexOf(next_gene, set, num_genes);
    exchanged        = set[to_exchange];

    // Construct the next solution
    if (gene_index != -1){
      set[gene_index]  = exchanged;
    }
    set[to_exchange] = next_gene;

    // Update the weights
    copyArr(weights, next_weights, t); // copy a -> b
    next_weights[set_to_exchange] = getWeight(set, k, set_to_exchange * k,
                                              n, A);
    other_set_to_exchange = gene_index / k;
    if (gene_index != -1 && set_to_exchange != other_set_to_exchange){
      next_weights[other_set_to_exchange] = getWeight(set, k,
                                            other_set_to_exchange * k, n, A);
    }
    next_weight = sum(next_weights, t);

    // Transition to the next state pseudorandomly based on the acceptance ratio
    log_acceptance_prob = c * next_weight - c * current_weight;
    logcoin = log((double) rand() / (double) RAND_MAX);
    if (logcoin <= log_acceptance_prob){
      if (gene_index != -1) num_within_swap += 1;
      current_weight = next_weight;
      copyArr(next_weights, weights, t);
      num_swapped += 1;
    }
    else{
      if (gene_index != -1){
        set[gene_index] = next_gene;
      }
      set[to_exchange] = exchanged; // Undo the swap
    }

    /* Thin out the chain by periodically saving solutions */
    if ( (i+1) % step_len == 0 ){
      for (j=0; j < num_genes; j++) gene_sets[step_index][j] = set[j];
      for (j=0; j < t; j++) set_weights[step_index][j] = weights[j];
      step_index += 1;
    }
    
    /* Output a simple progress bar */
    if ( (verbose == 1) && (i % prog_step == 0) ){
      steps = (int) 72. * i / num_iters;
      for (j = 0; j < steps; j++) printf("+");
      for (j = 0; j < 72 - steps; j++) printf(" ");
      printf(" [ %d%% ]\r", (int) (100. * i / num_iters));
      fflush(stdout);
    }
  }
  if (verbose == 1){
    printf("\n");
    printf("  - Swapped percentage: %d / %d (%d within)\n", num_swapped,
           num_iters, num_within_swap);
    printf("  - Look ups: %d\n", lookups);
    printf("  - Unique gene sets: %d\n", calculations);
    end = clock();
    printf("  - MCMC runtime: %f secs\n", (double)(end - begin) / CLOCKS_PER_SEC);
  }

  // Free memory
  free(set);
  free(weights);

}

/*****************************************************************************/
// Create a binary mutation matrix from the given mutation data
void init_mutation_data(bool **A, int m, int n, PyObject *gene2cases){
  // Declarations
  int i, j, num_cases, index;
  PyObject *cases, *val;

  // Iterate through every mutated sample in gene2cases to fill in A
  for (i = 0; i < m; i++){
    cases = PyList_GetItem(gene2cases, i);
    num_cases = (int) PyList_Size(cases);
    for (j = 0; j < num_cases; j++){
      val = PyList_GetItem(cases, j);
      index = (int) PyLong_AsLong(val);
      A[i][index] = true;
    }
  }

  // Free memory
  Py_DECREF(cases);
  Py_DECREF(val);

}

/*****************************************************************************/
// Functions callable from Python
long count_bits(long n) {     
  unsigned int c; // c accumulates the total bits set in v
  for (c = 0; n; c++) 
    n &= n - 1; // clear the least significant bit set
  return c;
}

PyObject *py_mcmc(PyObject *self, PyObject *args){
  // Parameters
  int i, j, ell;
  int k, t, step_len, num_iterations, numFrozen, numPatients, numGenes, verbose;
  bool **A;
  int **gene_sets;
  double **weights;
  PyObject *solns, *gene2cases, *soln, *val, *geneset;
  
  /* Parse Python arguments */
  if (! PyArg_ParseTuple( args, "iiiiO!iii", &k, &t, &numGenes, &numPatients,
                          &PyList_Type, &gene2cases, &num_iterations,
                          &step_len, &verbose)) {
    return NULL;
  }

  // Set up the mutation data as a giant mutation matrix
  A = malloc(sizeof(bool *) * numGenes);
  for (i = 0; i < numGenes; i++){
    A[i] = malloc(numPatients * sizeof(bool));
    for (j = 0; j < numPatients; j++) A[i][j] = false;
  }
  init_mutation_data(A, numGenes, numPatients, gene2cases);

  // Allocate memory for the gene sets
  numFrozen = num_iterations / step_len;
  gene_sets = malloc(sizeof(int *) * numFrozen);
  weights = malloc(sizeof(double *) * numFrozen);
  for (i=0; i < numFrozen; i++){
    gene_sets[i] = malloc(sizeof(int) * k * t);
    weights[i] = malloc(sizeof(double) * t);
  }
  
  // Run the MCMC
  mcmc(A, numGenes, numPatients, k, t, num_iterations, step_len,
       gene_sets, weights, verbose);
  
  // Convert the gene sets identified into Python objects
  solns = PyList_New(numFrozen);
  for (i = 0; i < numFrozen; i++){
    soln = PyList_New(t);
    for (j = 0; j < t; j++ ){
      geneset = PyList_New(k+1);
      for (ell = 0; ell < k; ell++){
        val = PyInt_FromLong((long) gene_sets[i][j * k + ell]);
        PyList_SET_ITEM(geneset, ell, val);
      }
      val = PyFloat_FromDouble(weights[i][j]);
      PyList_SET_ITEM(geneset, k, val);
      PyList_SET_ITEM(soln, j, geneset);
    }
    PyList_SET_ITEM(solns, i, soln);
  }

  // Free memory
  for (i = 0; i < numGenes; i++) free(A[i]);
  free(A);

  for (i = 0; i < numFrozen; i++){
    free(gene_sets[i]);
    free(weights[i]);
  }
  
  // Return results!
  return Py_BuildValue("O", solns); 
}

PyObject *py_exhaustive(PyObject *self, PyObject *args){
	return Py_BuildValue("");
}


/*****************************************************************************/
// Initialize the module
PyMethodDef MultiDendrixMethods[] = {
    {"mcmc", py_mcmc, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC initmultidendrixmodule(void) {
    PyObject *m = Py_InitModule("multidendrixmodule", MultiDendrixMethods);
    if (m == NULL) {
        return;
    }
}