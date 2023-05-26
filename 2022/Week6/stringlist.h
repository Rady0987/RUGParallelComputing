#pragma once

/**
 * This is a special list of strings implementation. It is quite tailored
 * to the usecases for this assignment. It will store a list of <count>
 * strings, each of maximum length <size>.
 * 
 * The stringList will have been allocated a single large block of contiguous
 * memory of <count> * <size> bytes. Then, every <size> bytes, a new string
 * can live. The <strings> list then contains pointers to each individual string.
 * 
 * In addition to a simple implementation of this list, we have provided some
 * functions that perform manipulations on all elements in this list, such as
 * making substitutions, uppercasing random letters in the list, etc. See below
 * for a full overview.
 */
struct stringList {
    char **strings;
    char *block;
    int count;
    int size;
};

struct stringList *allocStringList(int count, int stringSize);
void freeStringList(struct stringList *);

/**
 * This function is the most generic manipulation function we provide. It
 * allows one to substitute the <n>th occurence of the single character
 * <from> with string <to>, repeated for every string in the list.
 * 
 * It will skip empty strings and will not grow string beyond its maximum
 * size (from <size> in the stringList struct).
 */
struct stringList *manipulateList(struct stringList *, char from, char *to, int n);

/**
 * This function transforms every string in the list to its uppercase
 * variant. It will skip emtpy strings, but will not skip strings that
 * are the same in lower or uppercase (such as "123").
 */
struct stringList *uppercaseList(struct stringList *);

/**
 * This function generates a very long new list of length <size> * <count>
 * where each string in the original list is replicated <size> times, each
 * time with a different character capitalized.
 * 
 * It will skip empty strings and will generate extra empty strings when
 * the length of a particular string is less than <size>
 * 
 * For example: the string list "abc", "de" of <size> 3 and <count> 2 will
 * become "Abc", "aBc", "abC", "Def", "dEf", "" of <size> 3 and <count> 6,
 * where the last string is empty since "de" is shorter than <size> (3).
 * 
 */
struct stringList *capitalList(struct stringList *);

/**
 * This function generates a very long new list containing all possible
 * combinations of elements in the first and second list. The new list
 * will have size MAX(<size first>, <size second>), this allows you to
 * set a maximum password size on the source list and prevent any combination
 * that is longer than that.
 * 
 * The new list <count> will be <first count> * <second count>. Any empty
 * strings in the first list will be skipped (and leave <second count> emtpy
 * strings), similar for any empty strings in the second list.
 * 
 * For example: the string lists "abc", "def" and "zyx", "wvu", both of size
 * 6 and length 2, will generate the string list
 * "abczyx", "abcwvu", "defxyz", "defwvu" of size 6 and length 4.
 */
struct stringList *combinationList(struct stringList *, struct stringList *);

/**
 * This function will split a given stringlist into <n> equal parts. If
 * <count> is not a multiple of <n>, the last list will contain the
 * "leftover" items after making a nice division.
 */
struct stringList **splitList(struct stringList *, int n);
