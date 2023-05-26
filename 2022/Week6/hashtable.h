#pragma once

#include <stdlib.h>

#include "linkedlist.h"

/**
 * This is a string hashtable implementation. It is used to quickly deduplicate
 * lists of strings, but also to quickly look up any computed hashes to see if
 * they belong to a real user.
 * 
 * The hashtable must have a capacity that is a power of 2. (1024, 2048, etc)
 * (the allocate method will round up any capacity you hand it to the next
 * power of two - so you don't need to worry about this limitation).
 * For optimal performance, the hashtable has around double the capacity of the
 * expected number of unique keys it will store.
 * 
 * The structure is as follows: the table has a list of <capacity> entries. Each
 * entry is identified by a key (whose hash is the index in the array of entries)
 * and has a linked list as value - this allows one key to have multiple values.
 * Each item in this linked list is just a void-pointer, so it can point to anything
 * you like.
 */
struct hashTableEntry {
    char *key;
    struct listEntry *items;
};

struct hashTable {
    struct hashTableEntry *entries;
    size_t size;
    size_t capacity;
    size_t mask;
};

/**
 * Allocate a new hashtable with the given capacity
 */
struct hashTable *hashTableAlloc(size_t);

/**
 * Free a hashtable
 */
void hashTableFree(struct hashTable *);

/**
 * Insert a new item into the hashtable. The hash value of the key will be
 * determined, after which a new item will be stored in the linked list
 * associated with this key that will contain the void pointer you pass in here
 */
void hashTableInsert(struct hashTable *, char *, void *);

/**
 * Find the given key in the hashtable. Will return NULL if the key was
 * not found in the hashtable. Will return a pointer to the relevant
 * hashtTableEntry.
 */
struct hashTableEntry *hashTableGet(struct hashTable *, char *);
