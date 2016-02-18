#ifndef TESTS_H
#define TESTS_H

#include <iostream>
#include <string.h>
#include <iomanip>
#include "misc.h"

using namespace std;

#define POINTS 1024
#define PI atan(1) * 4

void fundamentalTest(unsigned int n);			// 07/11/2014
void test0();					// 27/06/2014
void test1();					// 22/07/2014
void test2(unsigned int n);		// 24/07/2014
int test3(const char *filepath, double percentage); // Compression test.
int test4(const char *filepath, double percentage); // Compression test.
int test5(float percentage);	// Intervalar and ponctual compression comparison.
int test6(float percentage);	// Matrix version of test5().
int test7(const char *filepath);	// Image quality test for ERAD 2015.
int test8(int n, int levels);	// Quality of the new Haar A-Trous implementarion.
int test9(int n, int max, int levels); // Quality of new Haar A-Trous 2D results.
int test10(const char *filepath, int type, int levels);	// Quality test for the new Haar A-Trous 2D (for papers).
int test11(const char *filepath, int type);	// Quality test for the new Haar Cascade 2D (for papers).

void test12Desc();				// Prints the description of test12, which is considered test 1.
void test12Param();				// Prints all parameters needed to perform test12.
int test12(int argc, char **argv);	// 1D cascade execution performance, based on test2.

void test13Desc();				// Prints the description of test13, which is considered test 2.
void test13Param();				// Prints all parameters needed to perform test13.
int test13(int argc, char **argv);	// 2D cascade execution performance, based on test12.

void test14Desc();				// Prints the description of test14, which is considered test 3.
void test14Param();				// Prints all parameters needed to perform test14.
int test14(int argc, char **argv);	// 1D cascade execution performance, based on test12.

void test15Desc();				// Prints the description of test15, which is considered test 4.
void test15Param();				// Prints all parameters needed to perform test15.
int test15(int argc, char **argv);	// 2D cascade execution performance, based on test13.

void matrixCopy(double **m1, double **m2, int x, int y);

typedef struct NewArgs
{
	int argc;
	char **argv;
} NewArgs;

NewArgs genNewArgs(int argc, char **argv, int takeOut);
void freeNewArgs(NewArgs& newArgs);
void listAllTests();
int testIndexer(int argc, char **argv);

#endif /* TESTS_H */