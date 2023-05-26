all: binaryTree rbTree

binaryTree: execute.c binaryTree.c
	gcc -O2 -fopenmp execute.c binaryTree.c -o binaryTree -lm

rbTree: execute.c rbTree.c
	gcc -O2 -fopenmp execute.c rbTree.c -o rbTree -lm
