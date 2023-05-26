#include <stdlib.h>

#include "linkedlist.h"

struct listEntry *listInit(void *value) {
    struct listEntry *list = malloc(sizeof(struct listEntry));
    list->next = list;
    list->prev = list;
    list->value = value;

    return list;
}

struct listEntry *listPush(struct listEntry *list, void *value) {
    struct listEntry *new = listInit(value);

    new->next = list;
    list->prev = new;

    return new;
}

struct listEntry *listPushBack(struct listEntry *list, void *value) {
    struct listEntry *new = listInit(value);

    new->prev = list;
    list->next = new;

    return new;
}

void *listPop(struct listEntry **list) {
    void *result = (*list)->value;

    if(listIsEnd(*list)) {
        free(*list);
        *list = NULL;
        return result;
    }

    *list = (*list)->next;
    free((*list)->prev);
    (*list)->prev = *list;

    return result;
}

void listFree(struct listEntry *list) {
    while(!listIsEnd(list)) {
        list = list->next;
        free(list->prev);
    }

    free(list);
}

void listPushP(struct listEntry **list, void *value) {
    *list = listPush(*list, value);
}

void listPushBackP(struct listEntry **list, void *value) {
    *list = listPushBack(*list, value);
}
