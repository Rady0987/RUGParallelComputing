#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "hashtable.h"
#include "linkedlist.h"

// FNV-1a hashing algorithm
#define FNV_PRIME 0x100000001b3ULL
#define FNV_OFFSET 0xcbf29ce484222325ULL
static uint64_t hashTableHash(char *string) {
    uint64_t result = FNV_OFFSET;
    unsigned char *c = (unsigned char *) string;

    while(*c) { // Loop until end of \0-terminated string
        result ^= *c;
        result *= FNV_PRIME;
        c++;
    }

    return result;
}

// Hash table operations
struct hashTable *hashTableAlloc(size_t size) {
    struct hashTable *table = malloc(sizeof(struct hashTable));

    // Our size might not be a power of two, but we need the capacity
    // to be a power of two. To find the nearest (upward) power of two:
    size_t capacity = 2; size--;
    while (size >>= 1) capacity <<= 1;

    if(capacity & (capacity - 1)) {
        // only accept powers of 2
        free(table);
        return NULL;
    } 
    table->size = 0;
    table->capacity = capacity;
    table->mask = capacity - 1;
    table->entries = calloc(capacity, sizeof(struct hashTableEntry));

    return table;
}

void hashTableFree(struct hashTable *table) {
    for(size_t i = 0; i < table->capacity; i++) {
        if(table->entries[i].key == NULL) continue;
        listFree(table->entries[i].items);
    }

    free(table->entries);
    table->entries = NULL;
    table->size = 0;
    table->capacity = 0;
    table->mask = 0;

    free(table);
}

void hashTableExpand(struct hashTable *table) {
    struct hashTableEntry *oldEntries = table->entries;
    size_t oldCapacity = table->capacity;

    // doubles capacity
    table->capacity *= 2;
    table->mask = table->capacity - 1;
    table->entries = calloc(table->capacity, sizeof(struct hashTableEntry));

    // re-set all data
    for(size_t i = 0; i < oldCapacity; i++) {
        if(oldEntries[i].key == NULL) continue;
        size_t idx = hashTableHash(oldEntries[i].key) & table->mask;

        while(table->entries[idx].key != NULL) {
            idx++;
            idx &= table->mask;
        }

        table->entries[idx].key = oldEntries[i].key;
        table->entries[idx].items = oldEntries[i].items;
    }

    free(oldEntries);
}

void hashTableInsert(struct hashTable *table, char *key, void *value) {
    // Expand if beyond 50%
    if(table->size * 2 > table->capacity)
        hashTableExpand(table);

    size_t i = hashTableHash(key) & table->mask;

    while(table->entries[i].key != NULL) {
        // If we find duplicate
        if(strcmp(key, table->entries[i].key) == 0) {
            listPushP(&(table->entries[i].items), value);
            return;
        }

        i++;
        i &= table->mask;
    }

    // No existing entry, create new entry at current index
    table->entries[i].key = key;
    table->entries[i].items = listInit(value);
    table->size++;
}

struct hashTableEntry *hashTableGet(struct hashTable *table, char *key) {
    size_t i = hashTableHash(key) & table->mask;

    while(table->entries[i].key != NULL) {
        // If we find the item
        if(strcmp(key, table->entries[i].key) == 0) {
            return &table->entries[i];
        }

        i++;
        i &= table->mask;
    }

    return NULL;
}
