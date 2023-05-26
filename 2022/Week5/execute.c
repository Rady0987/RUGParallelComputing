#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <sys/resource.h>

#include "tree.h"

/**
 * This file will handle all the I/O for the trees program, it will also set the
 * number of OpenMP threads already to OMP so you should be all set after this.
 * This code calls parallelSum() in tree.h, which is the point where you should
 * start with your parallel implementation.
 * 
 * Do NOT change this file; also do NOT submit this file; Themis will use a custom
 * setup that performs some additional performance checks and will have different
 * input/output.
 */

#define MAX_VALUE 1000000

#define RANDOM 0
#define SEQUENTIAL 1

long *initialise_values(int N, int mode) {
    long *arr = (long*) malloc(N*sizeof(long));
    
    if (mode == RANDOM) {
        for (int i = 0; i < N; i++) {
            arr[i] = (long) (rand() % MAX_VALUE);
        }
    }
    
    if (mode == SEQUENTIAL) {
        for (int i = 0; i < N; i++) {
            arr[i] = (long) i;
        }
    }
    
    return arr;
}

void treePerformance(Tree inputTree, int threads) {    
    double begin = omp_get_wtime();
    
    long result = parallelSum(inputTree, threads);
        
    double end = omp_get_wtime();

    double exec_time = (double) (end - begin);
    fprintf(stderr, "Execution time: %f seconds\n", exec_time);
    fprintf(stderr, "Result for tree sum = %ld\n", result);
}

/**
 * Code used for initialising the tree algorithms. Needs one argument: the number of threads.
 * Usage: ./a.out <num-threads>
 * The input will be read from stdin (depth and mode); depth is the tree depth and mode is 0
 * for RANDOM or 1 for SEQUENTIAL.
 */
int main(int argc, char *argv[]) {  
    if(argc != 2) {
        fprintf(stderr, "Fatal error: please provide the number of threads to run.\n");
        exit(EXIT_FAILURE);
    }

    int N, mode;
    int num_threads = atoi(argv[1]);
    scanf("%d %d", &N, &mode);

    // Set custom stack size
    struct rlimit rl;
    if((getrlimit(RLIMIT_STACK, &rl)) != 0) {
        fprintf(stderr, "Failure increasing stack size\n");
        exit(EXIT_FAILURE);
    }

    long desiredStackSize = 512000000l;
    rl.rlim_cur = (rl.rlim_max != RLIM_INFINITY && desiredStackSize > rl.rlim_max) ? rl.rlim_max : desiredStackSize;

    if(setrlimit(RLIMIT_STACK, &rl) != 0) {
        fprintf(stderr, "Failure increasing stack size\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Increased stack size to %ld.\n", rl.rlim_cur);

    // We seed the RNG here to make sure the output is reproducible on Themis
    srand(42);
    
    // Initialise number of threads for OMP
    omp_set_num_threads(num_threads);
    
    // Initialise values for both trees
    fprintf(stderr, "Initialising values...");
    long *vals = initialise_values(N, mode);
    fprintf(stderr, " DONE\n");
    
    // Compute true sum
    long sum = 0;
    for (int i = 0; i < N; i++) {
        sum += vals[i];
    }
    fprintf(stderr, "True sum = %ld\n", sum);
    
    // Construct the tree
    fprintf(stderr, "Constructing Tree...");
    Tree tree = (mode == SEQUENTIAL) ? constructSequentialTree(N) : constructFromArray(vals, N);
    int treeDepth = computeDepth(tree);
    fprintf(stderr, "tree depth = %d\n", treeDepth);
    
    // Analyse performance of the tree
    fprintf(stderr, "Evaluating tree summation...\n");
    treePerformance(tree, num_threads);
    
    // Free data
    free(vals);
    freeTree(tree);
    
    return 0;
}
