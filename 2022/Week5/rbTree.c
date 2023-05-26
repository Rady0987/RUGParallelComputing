#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <omp.h>

#include "tree.h"

/**
 * Implementation of an RB tree
 */

#define RIGHT 0
#define LEFT 1

#define BLACK 0
#define RED 1

#define SENTINEL_VALUE LONG_MIN 

// This Node implementation is a bit more complex than the one for the normal
// binary tree; it includes the parent and the colour.
struct Node {
    long value;
    int colour;
    Tree left, right;
    Tree parent;
};

Tree NIL = NULL;

// We predefine some functions here that are not included in the header file,
// but that are used by our implementation
Tree emptyTree();
Tree newTree(long val, Tree treeLeft, Tree treeRight);
void registerValue(Tree *root, long val);
long computeSumSerial(Tree root);


/**
 * Construction operations
 */

Tree createSentinel() {
    Tree sentinel = malloc(sizeof(struct Node));
    sentinel->value = SENTINEL_VALUE;
    sentinel->left = NULL;
    sentinel->right = NULL;
    sentinel->colour = BLACK;
    
    return sentinel;
}

void initialiseSentinel() {
    NIL = createSentinel();
}

Tree emptyTree() {
    if (NIL == NULL) {
        initialiseSentinel();
    }
    
    return NIL;
}

Tree newTree(long val, Tree treeLeft, Tree treeRight) {
    if (NIL == NULL) {
        initialiseSentinel();
    }
    
    Tree new = malloc(sizeof(*new));
    new->value = val;
    new->left = treeLeft;
    new->right = treeRight;
    new->parent = NULL;
    
    return new;
}

Tree constructFromArray(long *arr, int length) {
    if (NIL == NULL) {
        initialiseSentinel();
    }
    
    Tree new = malloc(sizeof(*new));
    new->value = arr[0];
    new->left = NIL;
    new->right = NIL;
    new->parent = NULL;
    new->colour = RED;
        
    for (int i = 1; i < length; i++) {
        registerValue(&new, arr[i]);
    }
    
    return new;
}

Tree constructSequentialTree(int N) {
    if (NIL == NULL) {
        initialiseSentinel();
    }

    Tree new = malloc(sizeof(*new));
    new->value = 0;
    new->left = NIL;
    new->right = NIL;
    new->parent = NULL;
    new->colour = RED;

    for (int i = 1; i < N; i++) {
        registerValue(&new, (long) i);
    }

    return new;
}

void freeTree(Tree root) {
    if (root->left != NIL) {
        freeTree(root->left);
    }
    if (root->right != NIL) {
        freeTree(root->right);
    }
    
    free(root);
}


/**
 *  Operations on tree
 */
void leftRotate(Tree root, Tree node) {
    Tree curr_node = node->right;
    
    // Turn the left sub-tree of curr_node into the right sub-tree of node
    node->right = curr_node->left;
    if (curr_node->left != NIL) {
        curr_node->left->parent = node;
    }
    
    curr_node->parent = node->parent;
    
    // Let the parent of node point to curr_node
    if (node->parent == NULL) {
        //printf("changing root\n");
        root = curr_node;
        curr_node->parent = NULL;
    }
    else {
        if (node == node->parent->left) {
            node->parent->left = curr_node;
        }
        else { // Node was on the right
            node->parent->right = curr_node;
        }
    }
    
    // Place node on curr_node's left
    curr_node->left = node;
    node->parent = curr_node;
}

void rightRotate(Tree root, Tree node) {
    Tree curr_node = node->left;
    
    // Turn the right sub-tree of curr_node into the left sub-tree of node
    node->left = curr_node->right;
    if (curr_node->right != NIL) {
        curr_node->right->parent = node;
    }
    
    curr_node->parent = node->parent;
    
    // Let the parent of node point to curr_node
    if (node->parent == NULL) {
        root = curr_node;
    }
    else {
        if (node == node->parent->right) {
            node->parent->right = curr_node;
        }
        else { // Node was on the left
            node->parent->left = curr_node;
        }
    }
    
    // Place node on curr_node's left
    curr_node->right = node;
    node->parent = curr_node;
}

void restore(Tree *root, Tree node) {    
    Tree curr_node;
    
    if (node->parent->parent == NULL) { // Special case: depth = 1
        node->colour = RED;
        (*root)->colour = BLACK;
        return;
    }
    
    while (node != *root && node->parent->colour == RED) {
        if (node->parent == node->parent->parent->left) {
            curr_node = node->parent->parent->right; // Right uncle
            
            if (curr_node->colour == RED) { // change to BLACK
                node->parent->colour = BLACK;
                curr_node->colour = BLACK;
                node->parent->parent->colour = RED;
                node = node->parent->parent; // move node up in tree
            }
            else { // curr_node is BLACK
                if (node == node->parent->right) {
                    node = node->parent;
                    leftRotate(*root, node);
                }
                node->parent->colour = BLACK;
                node->parent->parent->colour = RED;
                rightRotate(*root, node->parent->parent);
            }
        }
        else { // Same procedure, but switch left and right
            curr_node = node->parent->parent->left; // Left uncle
            
            if (curr_node->colour == RED) { // change to BLACK
                node->parent->colour = BLACK;
                curr_node->colour = BLACK;
                node->parent->parent->colour = RED;
                node = node->parent->parent; // move node up in tree
            }
            else { // curr_node is BLACK
                if (node == node->parent->left) {
                    node = node->parent;
                    rightRotate(*root, node);
                }
                node->parent->colour = BLACK;
                node->parent->parent->colour = RED;
                leftRotate(*root, node->parent->parent);
            }
        }
    }
    
    // Restore root to top
    while ((*root)->parent != NULL) {
        *root = (*root)->parent;
    }
    
    (*root)->colour = BLACK;
}

void registerValue(Tree *root, long val) {
    assert(root != NULL);
    
    Tree curr_node = *root;
    Tree prev_node = curr_node;
    int direction = -1;
        
    while (curr_node != NIL) {
        prev_node = curr_node;
        if (curr_node->value < val) {
            // Add to right side of sub-tree
            direction = RIGHT;
            curr_node = curr_node->right;
        } else {
            // Add to left side of sub-tree
            direction = LEFT;
            curr_node = curr_node->left;
        }
    }
        
    Tree new_node = newTree(val, NIL, NIL);
    
    // Add new node to the original tree
    if (direction == LEFT) {
        prev_node->left = new_node;
    } else {
        prev_node->right = new_node;
    }
    
    // Register parent
    new_node->parent = prev_node;
    
    // Assign colour RED
    new_node->colour = RED;
            
    // Restore red-black property
    restore(root, new_node);
}

int computeDepth(Tree root) {
    if (root == NIL) {
        return 0;
    }
    
    int depthLeft = computeDepth(root->left);
    int depthRight = computeDepth(root->right);
        
    return 1 + ( depthLeft > depthRight ? depthLeft : depthRight );
}

long computeSumSerial(Tree root) {
    if (root == NIL) {
        return 0;
    }
    
    long val = root->value;
    val += computeSumSerial(root->left);
    val += computeSumSerial(root->right);
    
    return val;
}

void printTree(Tree root) {
    if (root == NIL) {
        return;
    }
    
    printTree(root->left);
    printf("%ld\n", root->value);
    printTree(root->right);
}

long calculateSum(Tree root, int threads) {
    long result;

    if (root == NULL) {
        return 0;
    }

    result = root->value;

    /* The tasks are created at each depth level of the tree as long as the number of 
    threads is greater than one, at each recursive call, the number of available threads 
    is reduced by 1 */
    if(threads > 1) {

        /* Each task sets the result variable as shared, otherwise each 
        thread would create its own local variable, then throw away the results. */
        #pragma omp task shared(result)  
        result += calculateSum(root->left, threads - 1);

        #pragma omp task shared(result)
        result += calculateSum(root->right, threads - 1);
    
        /* This barrier ensures that the program flow will get paused until all queued 
        tasks have been executed. */
        #pragma omp taskwait    

    } else {
        /* At the moment when there is only one available thread, the sum of the remaining 
        tree nodes are calculated serially */
        return computeSumSerial(root);
    }

    return result;
}

long parallelSum(Tree root, int threads) {
    long result;

    /* The parallel pragma clause all of the threads in the pool to execute the next 
    block of code. The single directive causes all threads but one (usually the first to encounter 
    the block) to not execute it, while the nowait turns off the barrier on the single; there's already 
    a barrier on the enclosing parallel region, to which the other threads will rush. */
    #pragma omp parallel    
    #pragma omp single nowait
    result = calculateSum(root, threads);

    return result;
}
