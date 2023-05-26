#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "linkedlist.h"
#include "hashtable.h"
#include "stringmanipulation.h"
#include "stringlist.h"
#include "pw_helpers.h"
#include <mpi.h>

/**
 * Our entrypoint. We require two arguments to our program: the paths to a passwd and
 * shadow file. The number of threads/processes is dictated by MPI, and is out of our
 * control at this point.
 * 
 * Run like: mpiexec -n <threads> ./guessword <passwd> <shadow>
 */
int main(int argc, char **argv) {
    // Check arguments
    if(argc != 3) {
        fprintf(stderr, "Usage: ./guessword <passwd> <shadow>");
        return EXIT_FAILURE;
    }

    MPI_Init(&argc, &argv);

    int numProcesses, processId;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    ///////////////////////////////////////////////////////////////////
    // We now set up the local environment
    ///////////////////////////////////////////////////////////////////
    
    // Read the password/shadow files and parse all input
    char *passwdPath = argv[1];
    char *shadowPath = argv[2];

    struct users users = parseInput(passwdPath, shadowPath, false);

    // Read precomputed guess list
    struct stringList *pwListMain = readStringsFile("Files/top250.txt", MAX_PW_LENGTH);
    tryPasswords(pwListMain, users.passwords, users.hashSetting);

    struct stringList *appendedTop250 = manipulateList(users.usernames, 'o', "[]", 1);
    tryPasswords(appendedTop250, users.passwords, users.hashSetting);
    free(appendedTop250);

    ///////////////////////////////////////////////////////////////////
    // We will now start to do the real work
    ///////////////////////////////////////////////////////////////////

    // We will try the provided list of passwords and all usernames appended
    // with 00.
    struct stringList *appendedPasswords = manipulateList(users.usernames, '\0', "00", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);
    free(appendedPasswords);

    // We will try the provided list of passwords and all first names with []
    // instead of o.
    appendedPasswords = manipulateList(users.names[0], 'o', "[]", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);
    free(appendedPasswords);

    // We will try the provided list of passwords and all first names with (
    // instead of c.
    appendedPasswords = manipulateList(users.names[0], 'c', "(", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);
    free(appendedPasswords);

    // We will try the provided list of passwords and all first names with [\]
    // instead of n.
    appendedPasswords = manipulateList(users.names[0], 'n', "[\\]", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);
    free(appendedPasswords);

    // We will try the provided list of passwords and all first names with [)
    // instead of d.
    appendedPasswords = manipulateList(users.names[0], 'd', "[)", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);
    free(appendedPasswords);

    // We will try the provided list of passwords and all first names with &
    // instead of e.
    appendedPasswords = manipulateList(users.names[0], 'e', "&", 1);
    tryPasswords(appendedPasswords, users.passwords, users.hashSetting);

    // We will try the first name, last name and middle name and see if they
    // match with the passwords.
    tryPasswords(users.names[1], users.passwords, users.hashSetting);
    tryPasswords(users.names[2], users.passwords, users.hashSetting);
    tryPasswords(users.names[0], users.passwords, users.hashSetting);


    ///////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////

    // Clean password list
    freeStringList(pwListMain);
    freeStringList(appendedPasswords);

    // Free users struct/information
    freeUserData(users);

    MPI_Finalize();
}