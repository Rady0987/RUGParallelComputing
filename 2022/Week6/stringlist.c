#include <stdlib.h>

#include "stringlist.h"
#include "stringmanipulation.h"

struct stringList *allocStringList(int count, int stringSize) {
    struct stringList *result = malloc(sizeof(*result));
    char **list = calloc(count, sizeof(*list));
    char *block = calloc(stringSize * count, sizeof(*block));

    for(int i = 0; i < count; i++) {
        list[i] = &block[stringSize * i];
    }

    *result = (struct stringList){
        .strings = list,
        .block = block,
        .count = count,
        .size = stringSize
    };
    return result;
}

void freeStringList(struct stringList *list) {
    free(list->block);
    free(list->strings);
    free(list);
}

struct stringList *manipulateList(struct stringList *list, char from, char *to, int n) {
    struct stringList *newList = allocStringList(list->count, list->size);

    for(int i = 0; i < list->count; i++) {
        if(strempty(list->strings[i])) continue; // skip

        manipulate(list->strings[i], from, to, &newList->strings[i], list->size, n);
    }

    return newList;
}

struct stringList *uppercaseList(struct stringList *list) {
    struct stringList *newList = allocStringList(list->count, list->size);

    for(int i = 0; i < list->count; i++) {
        if(strempty(list->strings[i])) continue; // skip

        strcpy(newList->strings[i], list->strings[i]);
        strtoupper(newList->strings[i]);
    }

    return newList;
}

struct stringList *capitalList(struct stringList *list) {
    struct stringList *newList = allocStringList(list->count * list->size, list->size);

    for(int i = 0; i < list->count; i++) {
        if(strempty(list->strings[i])) continue; // skip

        int len = strlen(list->strings[i]);

        for(int j = 0; j < len; j++) {
            strcpy(newList->strings[i * list->size + j], list->strings[i]);
            chartoupper(&newList->strings[i * list->size + j][j]);
        }
    }

    return newList;
}

struct stringList *combinationList(struct stringList *first, struct stringList *second) {
    int maxSize = first->size > second->size ? first->size : second->size;
    struct stringList *list = allocStringList(first->count * second->count, maxSize);

    for(int i = 0; i < first->count; i++) {
        if(strempty(first->strings[i])) continue; // skip

        int len = strlen(first->strings[i]);

        for (int j = 0; j < second->count; j++) {
            if(strempty(second->strings[j])) continue; // skip

            strcpy(list->strings[i * second->count + j], first->strings[i]);
            strcpy(list->strings[i * second->count + j] + len, second->strings[j]);
        }
    }

    return list;
}

struct stringList **splitList(struct stringList *list, int parts) {
    struct stringList **result = malloc(parts * sizeof(*result));
    int size = list->count / parts;
    int leftover = list->count % parts;

    for(int i = 0; i < parts; i++) {
        if(i == parts - 1) {
            result[i] = allocStringList(size + leftover, list->size);
        } else {
            result[i] = allocStringList(size, list->size);
        }

        memcpy(result[i]->block, list->strings[i * size], result[i]->count * result[i]->size);
    }

    return result;
}
