#pragma once

/**
 * This is a generic header file that is shared for both the normal binary
 * tree and the RB tree. The header file is implementation-agnostic and does
 * not care about how our "struct Node" is implemented; it just knows we have
 * a pointer to it. This allows us to swap out the implementation when compiling
 * and saves us from having to write two sets of code.
 */

typedef struct Node *Tree;

// Tree construction
Tree constructFromArray(long *arr, int length);
Tree constructSequentialTree(int N);
void freeTree(Tree root);

// Operations on tree
int computeDepth(Tree root);

// Parallel sum
// This function will be called by Themis/execute.c, this is where the parallel evaluation
// should start.
long parallelSum(Tree root, int threads);
