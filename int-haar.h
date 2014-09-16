#include <iostream>
#include <cmath>
#include <interval.hpp>

using namespace std;
using namespace cxsc;

#ifndef UINT
#define UINT unsigned int
#endif

#define HAAR_COMPRESS_ERROR 0.0000000001

#define INTHAAR
#define HAAROPMIZATION

void Haar_Composition(double *vec, int n, bool normal);
void Haar_CompositionStep(double *vec, int n, bool normal);
void Haar_Compress(double *vec, int n, float percentage);
void Haar_Levels_Compress(double *vec, int n, double percentage);
void Haar_Decomposition(double *vec, int n, bool normal);
void Haar_DecompositionStep(double *vec, int n, bool normal);
double **Haar_Decomposition_For_Graphs(double *vec, int n, bool normal);
void Haar_MatrixComposition(double **matrix, int rows, int cols, bool normal, bool standard);
void Haar_MatrixDecomposition(double **matrix, int rows, int cols, bool normal, bool standard);
void Haar_NonStandardComposition(double **matrix, int rows, int cols, bool normal);
void Haar_NonStandardDecomposition(double **matrix, int rows, int cols, bool normal);
void Haar_StandardComposition(double **matrix, int rows, int cols, bool normal);
void Haar_StandardDecomposition(double **matrix, int rows, int cols, bool normal);

#ifdef HAAROPMIZATION

void VinisMatrixNormalization(double **mat, UINT n, bool standard, bool invert = false);
void VinisNonStandardMatrixNormalization(double **matrix, UINT n, bool invert = false);
void VinisNormalization(double *vec, UINT n);
void VinisStandardMatrixNormalization(double **mat, UINT n, bool invert = false);

#endif /* HAAROPMIZATION */

#ifdef INTHAAR

void INT_Haar_Composition(interval *vec, int n, bool normal);
void INT_Haar_CompositionStep(interval *vec, int n, bool normal);
void INT_Haar_Decomposition(interval *vec, int n, bool normal);
void INT_Haar_DecompositionStep(interval *vec, int n, bool normal);
void INT_Haar_MatrixComposition(interval **matrix, int rows, int cols, bool normal, bool standard);
void INT_Haar_MatrixDecomposition(interval **matrix, int rows, int cols, bool normal, bool standard);
void INT_Haar_NonStandardComposition(interval **matrix, int rows, int cols, bool normal);
void INT_Haar_NonStandardDecomposition(interval **matrix, int rows, int cols, bool normal);
void INT_Haar_StandardComposition(interval **matrix, int rows, int cols, bool normal);
void INT_Haar_StandardDecomposition(interval **matrix, int rows, int cols, bool normal);

#ifdef HAAROPMIZATION

void INT_VinisMatrixNormalization(interval **mat, UINT n, bool standard, bool invert = false);
void INT_VinisNonStandardMatrixNormalization(interval **matrix, UINT n, bool invert = false);
void INT_VinisNormalization(interval *vec, UINT n);
void INT_VinisStandardMatrixNormalization(interval **mat, UINT n, bool invert = false);

#endif /* HAAROPMIZATION */

real INT_diameter(interval x);
real INT_error(interval *x, int n);
real INT_error(interval **x, int linhas, int colunas);

#endif /* INTHAAR */
