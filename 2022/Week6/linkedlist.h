#pragma once

#include <stdbool.h>

/**
 * This is a very simple linked list implementation. It is quite generic since
 * its value is just some void pointer, and thus can point to anything (just
 * keep track of the real type yourself and cast the pointer to the right type).
 * 
 * The convention in this implementation is that for the first item, the prev-
 * pointer will point to itself, and similarly for the last item, the next-pointer
 * will point to itself. This allows an easy implementation of listIsStart and
 * listIsEnd.
 */
struct listEntry {
    struct listEntry *prev;
    struct listEntry *next;
    void *value;
};

static inline bool listIsStart(struct listEntry *list) {
    return list == list->prev;
}

static inline bool listIsEnd(struct listEntry *list) {
    return list == list->next;
}

struct listEntry *listInit(void *);

// Add a new item to the beginning of the linked list, where the given listEntry
// is assumed to be the first entry of the list.
struct listEntry *listPush(struct listEntry *, void *);

// Add a new item to the end of the linked list, where the given listEntry is
// assumed to be the LAST entry of the list. This function will not check this!
struct listEntry *listPushBack(struct listEntry *, void *);

// Return the first item from the front of the list. If this is the only item
// in the list, the list will be freed.
void *listPop(struct listEntry **);

// Free the entire list
void listFree(struct listEntry *);

// Versions of listPush and listPushBack that operate in-place; you pass a pointer
// to the list pointer and it will modify the list in-place - this way you do not
// have to fiddle with return values.
void listPushP(struct listEntry **, void *);
void listPushBackP(struct listEntry **, void *);
