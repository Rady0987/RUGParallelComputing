#pragma once

#include <ctype.h>
#include <string.h>
#include <stddef.h>
#include <stdbool.h>

/**
 * This file contains some useful/small string manipulation operations.
 * These are mainly used by the stringList manipulations in stringlist.c.
 * 
 * It also contains some useful helper functions.
 */

// Determine whether a string is empty (pointer to NULL or first character
// being EOF \0 directly).
static inline bool strempty(char *s) {
    return s == NULL || *s == '\0';
}

// Perform a safe version of strdup() that "passes through" any NULLs.
static inline char *strsafedup(char *s) {
    return (s == NULL) ? NULL : strdup(s);
}

// Transforms a string to its lowercase variant in-place, so it modifies
// the original string. In addition, it will return the pointer to the
// string itself to facilitate function call nesting.
static inline char *strtolower(char *s) {
    if(s == NULL) return NULL;
    char *copy = s;
    for(; *copy; copy++) *copy = tolower(*copy);
    return s;
}

// Transforms a string to its uppercase variant in-place, so it modifies
// the original string. In addition, it will return the pointer to the
// string itself to facilitate function call nesting.
static inline char *strtoupper(char *s) {
    if(s == NULL) return NULL;
    char *copy = s;
    for(; *copy; copy++) *copy = toupper(*copy);
    return s;
}

// Variant of toupper() that changes the given character pointer in-place,
// so you don't have to deal with the return value.
static inline void chartoupper(char *s) {
    *s = toupper(*s);
}

// Variant of tolower() that changes the given character pointer in-place,
// so you don't have to deal with the return value.
static inline void chartolower(char *s) {
    *s = tolower(*s);
}

/**
 * Main string manipulation function that allows one to substitute the <n>th
 * occurrence of the single character <from> with the string <to>. It will
 * store the new string in the memory pointed to by result. The function
 * will make sure to not generate strings longer than maxSize.
 * 
 * If a string manipulation would have resulted in a string that is too long,
 * no action will be taken and the result buffer will remain untouched.
 */
static inline void manipulate(char *s, char from, char *to, char **result, size_t maxSize, int n) {
    if(strempty(s) || strlen(s) + strlen(to) > maxSize) {
        return;
    }

    char *occurrence = s - 1;
    for(int i = 0; i < n; i++) {
        occurrence = strchr(occurrence + 1, from);

        if(occurrence == NULL) {
            *result = NULL;
            return;
        }
    }

    ptrdiff_t before_len = (occurrence - s);
    ptrdiff_t after_len = strlen(s) - before_len; // one too long to include terminating \0

    memcpy(*result, s, before_len);
    memcpy((*result) + before_len, to, strlen(to));
    memcpy((*result) + before_len + strlen(to), occurrence + 1, after_len);
}
