#pragma once

#include <stdlib.h>
#include <stdbool.h>
#include "stringlist.h"
#include "hashtable.h"

// This is our personal Parallel Programming guarantee
#define MAX_PW_LENGTH 32
#define MAX_HASH_LENGTH 30

/**
 * Struct representing the data of a single user; where
 * names[4] is an array of 4 pointers to strings, of which
 * some elements might be NULL. This contains the user full
 * name split on the space. So, "Arnold Meijster" would become
 * [ (pointer to) "Arnold", (pointer to) "Meijster", NULL, NULL ].
 */
struct userData {
    int id;
    char *username;
    char *names[4];
    char *passwordHash;
};

/**
 * Struct representing all the users in the input files. There
 * will be <count> entries in each array. The <users> field is
 * a list of <count> struct userData entries, representing the
 * data of each user individually. Then, the usernames field is
 * a "string list" containing all UNIQUE, LOWERCASE usernames.
 * The string list has <count> entries, but any duplicate usernames
 * are set to the empty string.
 * 
 * Similary, the <names> field is again a list of 4 string lists,
 * corresponding to the list of all UNIQUE, LOWERCASE first names,
 * second names, etc - just like in the struct userData above.
 * Lastly, the passwords field is a hashtable containing all password
 * hashes to guess. This is used by tryPasswords to quickly lookup
 * the computed hashes.
 */
struct users {
    int count;
    struct userData *users;
    struct stringList *names[4];
    struct stringList *usernames;

    char *hashSetting; // used for crypt_r()
    struct hashTable *passwords;
};

/**
 * This function will parse the contents of the given files. It will return
 * all information from these files in the form of the structs defined
 * above. Very convenient. If silent is true, this will not print debugging
 * output.
 */
struct users parseInput(char *passwdPath, char *shadowPath, bool silent);

// To free the allocated data in parseInput
void freeUserData(struct users);

/**
 * This function will read a file of strings to a struct stringList, which
 * can be used to easily and quickly read in additional data. The file should
 * simply contain all entries on separate lines, separated by just \n. Each
 * entry should not be longer than the provided maxLength parameter (which
 * you should keep track of yourself).
 */
struct stringList *readStringsFile(char *filename, size_t maxLength);

/**
 * This function will hash all the passwords in the given stringList, will
 * look them up in the provided passwords hashtable, and if a result is found
 * it will print the result in the correct output format to the console without
 * delay.
 */
void tryPasswords(struct stringList *guesses, struct hashTable *passwords, char *hashSetting);
