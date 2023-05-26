#define _DEFAULT_SOURCE
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pwd.h>
#include <shadow.h>
#include <errno.h>
#include <crypt.h>

#include "pw_helpers.h"
#include "stringlist.h"
#include "stringmanipulation.h"
#include "hashtable.h"

// We assume the maximum length of the information for any user is 1024 bytes.
// That is fine for this assignment.
#define USER_BUFFER_LENGTH 1024

// We use a helper function to determine the number of lines in a file
#define LINECOUNT_BUFFER_SIZE 65536
int lineCount(FILE *file) {
    // We read the file in batches
    char buffer[LINECOUNT_BUFFER_SIZE];
    int count = 0;

    // Set the read position to the start of the file
    fseek(file, 0, SEEK_SET);

    while(1) {
        size_t bytesRead = fread(buffer, 1, LINECOUNT_BUFFER_SIZE, file);
        if(bytesRead == 0) break;

        for(size_t i = 0; i < bytesRead; i++)
            if(buffer[i] == '\n') count++;
    }

    // Reset the read position again
    fseek(file, 0, SEEK_SET);

    return count;
}


/**
 * This function will parse the contents of the given files. It will return
 * all information from these files in the form of the structs defined
 * above. Very convenient. If silent is true, this will not print debugging
 * output.
 */
struct users parseInput(char *passwdPath, char *shadowPath, bool silent) {
    // Sanity check
    if(passwdPath == NULL || shadowPath == NULL) {
        if(!silent) fprintf(stderr, "WARNING: Input file names are empty!");
        return (struct users) { .count = 0, .users = NULL, .names = {NULL, NULL, NULL, NULL}, .usernames = NULL };
    }

    if(!silent) fprintf(stderr, "Reading password files... ");

    // Open files for reading
    FILE *shadow = fopen(shadowPath, "r");
    FILE *passwd = fopen(passwdPath, "r");

    if(passwd == NULL || shadow == NULL) {
        if(!silent) fprintf(stderr, "WARNING: Could not open input files!");
        return (struct users) { .count = 0, .users = NULL, .names = {NULL, NULL, NULL, NULL}, .usernames = NULL };
    }

    // Determine the file length
    int count = lineCount(shadow);

    // Allocate space for the result struct. We use a "compound literal" here,
    // to easily and quickly initialise a struct.
    struct users users = (struct users) {
        .count = count,
        .users = calloc(count, sizeof(struct userData)),
        .hashSetting = NULL,
        .usernames = allocStringList(count, MAX_PW_LENGTH),
        .passwords = hashTableAlloc(2 * count), // 2x too large for better performance
        .names = {
            [0] = allocStringList(count, MAX_PW_LENGTH),
            [1] = allocStringList(count, MAX_PW_LENGTH),
            [2] = allocStringList(count, MAX_PW_LENGTH),
            [3] = allocStringList(count, MAX_PW_LENGTH),
        }
    };

    ///////////////////////////////////////////////////////////////////
    // We will now read all user data from the files provided. We do so
    // using the appropriate libraries that perform all this work for
    // us, which is very convenient. We will parse and process this data
    // and put it in our custom structs.
    ///////////////////////////////////////////////////////////////////

    // Allocate some variables that will contain the information read from the file
    // by the system libraries.
    struct passwd userinfoEntryBuffer, *userinfoEntry;
    char userinfoStringBuffer[USER_BUFFER_LENGTH];

    // Read data for each user
    for(int i = 0; i < count; i++) {
        // Get user info entry
        int result = fgetpwent_r(passwd, &userinfoEntryBuffer, userinfoStringBuffer, USER_BUFFER_LENGTH, &userinfoEntry);
        if(result == ENOENT || result == ERANGE || userinfoEntry == NULL) break; // we assume this to be end of file

        struct spwd *pwdEntry = fgetspent(shadow);
        if(pwdEntry == NULL) break; // we assume this to be end of file as well

        // Fill user struct, again using a "compound literal"
        // We use strdup() to quickly clone a string without having to manually
        // allocate memory and copy over data. We need to do this because the
        // calls above (fgetpwent_r and fgetspent) will reuse the same memory
        // for the next user, which will overwrite the current data.
        users.users[i] = (struct userData) {
            .id = i,
            .username = strdup(userinfoEntry->pw_name),
            .passwordHash = strdup(pwdEntry->sp_pwdp),
            .names = {
                // We use "designated initializers" to set specific array values here
                [0] = strtolower(strsafedup(strtok(userinfoEntry->pw_gecos, " ,"))),
                [1] = strtolower(strsafedup(strtok(NULL, " ,"))),
                [2] = strtolower(strsafedup(strtok(NULL, " ,"))),
                [3] = strtolower(strsafedup(strtok(NULL, " ,")))
            }
        };

        // Finding the hashSetting information if that has not been done. This can
        // be done once for the entire file, since it is guaranteed that all entries
        // in a single file have the same salt. (it will be different between files!)
        if(users.hashSetting == NULL) {
            char *scheme = strtok(pwdEntry->sp_pwdp, "$"); // first item is the scheme
            char *salt = strtok(NULL, "$"); // second item is the salt

            asprintf(&users.hashSetting, "$%s$%s$", scheme, salt);
        }
    }

    if(!silent) fprintf(stderr, "Read %d users.\n", users.count);

    ///////////////////////////////////////////////////////////////////
    // We will now make the lists of all unique usernames, first names,
    // second names, etc. We do this by keeping track of a hashtable
    // of all names seen so far. If a name is already in this table,
    // it will not be put in the list of all unqiue names again.
    ///////////////////////////////////////////////////////////////////
    
    // We use one hashtable per names list; as we only need them to be unique within
    // the list - not across lists.
    struct hashTable *nameTables[4] = {
        [0] = hashTableAlloc(count),
        [1] = hashTableAlloc(count),
        [2] = hashTableAlloc(count),
        [3] = hashTableAlloc(count)
    };
    struct hashTable *usernameTable = hashTableAlloc(count);

    if(!silent) fprintf(stderr, "Deduplicating names... ");

    for(int i = 0; i < count; i++) {
        // Username; if it does not exist in the hash table it is a new
        // username so we add it to the list.
        if(users.users[i].username != NULL && 
                hashTableGet(usernameTable, users.users[i].username) == NULL) {
            hashTableInsert(usernameTable, users.users[i].username, NULL);
            strcpy(users.usernames->strings[i], users.users[i].username);
        }

        // We now do the same as above for all 4 names.
        for(int j = 0; j < 4; j++) {
            if(users.users[i].names[j] != NULL && 
                    hashTableGet(nameTables[j], users.users[i].names[j]) == NULL) {
                hashTableInsert(nameTables[j], users.users[i].names[j], NULL);
                strcpy(users.names[j]->strings[i], users.users[i].names[j]);
            }
        }
    }

    if(!silent) fprintf(stderr, "Done\n");

    ///////////////////////////////////////////////////////////////////
    // We will now build the hashtable of all password hashes (so...
    // a hashhashtable?).
    ///////////////////////////////////////////////////////////////////

    if(!silent) fprintf(stderr, "Building password hashtable... ");
    for(int i = 0; i < count; i++) {
        // Each entry in the hashtable contains a pointer to our userData struct for later reference
        hashTableInsert(users.passwords, users.users[i].passwordHash, &users.users[i]);
    }
    if(!silent) fprintf(stderr, "Inserted %d passwords\n", count);

    ///////////////////////////////////////////////////////////////////
    // Clean up
    ///////////////////////////////////////////////////////////////////
    
    hashTableFree(usernameTable);
    for(int j = 0; j < 4; j++) hashTableFree(nameTables[j]);

    // Close files
    fclose(shadow);
    fclose(passwd);

    return users;
}

// To free the allocated data in parseInput
void freeUserData(struct users users) {
    freeStringList(users.usernames);
    hashTableFree(users.passwords);
    for(int i = 0; i < 4; i++) freeStringList(users.names[i]);

    for(int i = 0; i < users.count; i++) {
        free(users.users[i].username);
        free(users.users[i].passwordHash);
        for(int j = 0; j < 4; j++) free(users.users[i].names[j]);
    }
    free(users.users);
    free(users.hashSetting);
}

/**
 * This function will read a file of strings to a struct stringList, which
 * can be used to easily and quickly read in additional data. The file should
 * simply contain all entries on separate lines, separated by just \n. Each
 * entry should not be longer than the provided maxLength parameter (which
 * you should keep track of yourself).
 */
struct stringList *readStringsFile(char *filename, size_t maxLength) {
    // Determine the number of entries
    FILE *file = fopen(filename, "r");
    if(file == NULL) {
        fprintf(stderr, "Couldn't open file %s!\n", filename);
        return allocStringList(1, maxLength);
    }

    int count = lineCount(file);

    // Put results in a stringList
    struct stringList *results = allocStringList(count, maxLength);

    for(int i = 0; i < count; i++) {
        fgets(results->strings[i], maxLength, file);
        strtok(results->strings[i], "\n"); // get rid of the newline at the end
    }

    fclose(file); // clean up
    return results;
}

/**
 * This function will hash all the passwords in the given stringList, will
 * look them up in the provided passwords hashtable, and if a result is found
 * it will print the result in the correct output format to the console without
 * delay.
 */
void tryPasswords(struct stringList *guesses, struct hashTable *passwords, char *hashSetting) {
    // Pre-allocate a crypt_data struct for our hashing results
    struct crypt_data cryptData;

    // Try all passwords in the list
    for(int i = 0; i < guesses->count; i++) {
        if(strempty(guesses->strings[i])) continue;

        // Hash the attempted password with the proper settings
        char *hash = crypt_r(guesses->strings[i], hashSetting, &cryptData);

        // Lookup the hash in the passwords list
        struct hashTableEntry *result;
        if((result = hashTableGet(passwords, hash)) == NULL) continue; // Not found

        // The users that have the password <guesses->strings[i]> are now in the
        // result->items linked list. Traverse that list and print all results.
        struct listEntry *currentUser;
        currentUser = result->items;
        while(1) {
            // Print the user information
            printf("%s:%s\n", ((struct userData *) currentUser->value)->username, guesses->strings[i]);
            fflush(stdout); // Make sure to immediately print and not buffer the results and have
                            // them be lost.

            if(listIsEnd(currentUser)) break;
            currentUser = currentUser->next;
        }
    }
}
