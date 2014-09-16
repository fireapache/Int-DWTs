#include <iostream>
#include <cmath>
#include <interval.hpp>
#include "common.h"

using namespace std;
using namespace cxsc;

#define HAAR_COMPRESS_ERROR 0.0000000001

#define INTHAAR
#define HAAROPMIZATION

void Haar_Composition(double *vec, int n, bool normal);
void Haar_CompositionStep(double *vec, int n, bool normal);
void Haar_Compression(double *vec, int n, float percentage);
void Haar_PerLevel_Compression(double *vec, int n, double percentage);
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

void VinisMatrixNormalization(double **mat, uint n, bool standard, bool invert = false);
void VinisNonStandardMatrixNormalization(double **matrix, uint n, bool invert = false);
void VinisNormalization(double *vec, uint n);
void VinisStandardMatrixNormalization(double **mat, uint n, bool invert = false);

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

void INT_VinisMatrixNormalization(interval **mat, uint n, bool standard, bool invert = false);
void INT_VinisNonStandardMatrixNormalization(interval **matrix, uint n, bool invert = false);
void INT_VinisNormalization(interval *vec, uint n);
void INT_VinisStandardMatrixNormalization(interval **mat, uint n, bool invert = false);

#endif /* HAAROPMIZATION */

real INT_diameter(interval x);
real INT_error(interval *x, int n);
real INT_error(interval **x, int linhas, int colunas);

#endif /* INTHAAR */
