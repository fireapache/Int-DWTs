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

void matrixCopy(double **m1, double **m2, int x, int y);

#endif /* TESTS_H */