#ifndef INTHAAR_H
#define INTHAAR_H

// TODO: Fix 1D intervalar compression!

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
void Haar_Matrix_Compression(double **matrix, int n, float percentage);
void Haar_PerLevel_Matrix_Compression(double **matrix, int n, float percentage);
void Haar_Level_Matrix_Compression(double **matrix, int n, float percentage);
void Haar_Decomposition(double *vec, int n, bool normal);
void Haar_DecompositionStep(double *vec, int n, bool normal);
double **Haar_Decomposition_For_Graphs(double *vec, int n, bool normal);
void Haar_MatrixComposition(double **matrix, int rows, int cols, bool normal, bool standard);
void Haar_MatrixDecomposition(double **matrix, int rows, int cols, bool normal, bool standard);
void Haar_NonStandardComposition(double **matrix, int rows, int cols, bool normal);
void Haar_NonStandardDecomposition(double **matrix, int rows, int cols, bool normal);
void Haar_StandardComposition(double **matrix, int rows, int cols, bool normal);
void Haar_StandardDecomposition(double **matrix, int rows, int cols, bool normal);

template <typename T>
T** Haar_atrous_Decomposition(T *vec, int n, int levels, bool normal)
{
	if (n <= 0 || levels <= 0) return NULL;

	T **data = new T*[levels * 2];
	T *D0, *D1;
	T divisor;

	for (int i = 0; i < levels * 2; ++i)
	{
		data[i] = new T[n];
	}

	// Setting normalization factor, if needed.
	if (normal) divisor = sqrt(T(2.0));
	else divisor = T(2.0);

	D0 = vec;

	// Calculating degree coefficients.
	for (int i = 0; i < levels; ++i)
	{
		for (int j = 0; j < n - 1; ++j)
		{
			data[i * 2][j] = (D0[j] + D0[j + 1]) / divisor;
		}

		data[i * 2][n - 1] = (D0[n - 1] + D0[0]) / divisor;

		D0 = data[i * 2];
	}

	D0 = vec;
	D1 = data[0];

	// Calculating wavelet coefficients.
	for (int i = 0; i < levels; ++i)
	{
		for (int j = 0; j < n ; ++j)
		{
			data[i * 2 + 1][j] = D0[j] - D1[j];
		}

		if (i + 1 < levels);
		{
			D0 = data[i * 2];
			D1 = data[(i + 1) * 2];
		}
	}

	return data;
}

template <typename T>
void Haar_atrous_Normalization(T *vec, T **data, int n, int levels, bool invert = false)
{
	if (n <= 0 || levels <= 0) return;

	T factor, *D0, *D1;

	// Normalizing degree coefficients.
	for (int i = 0; i < levels; ++i)
	{
		// Normalization factor rule.
		factor = pow(T(2.0), T(i + 1) / T(2.0));

		if (invert)
		{
			for (int j = 0; j < n; ++j)
			{
				data[i * 2][j] /= factor;
			}
		}
		else
		{
			for (int j = 0; j < n; ++j)
			{
				data[i * 2][j] *= factor;
			}
		}
	}

	D0 = vec;
	D1 = data[0];

	// Recalculating wavelet coefficients.
	for (int i = 0; i < levels; ++i)
	{
		for (int j = 0; j < n ; ++j)
		{
			data[i * 2 + 1][j] = D0[j] - D1[j];
		}

		if (i + 1 < levels);
		{
			D0 = data[i * 2];
			D1 = data[(i + 1) * 2];
		}
	}
}

#ifdef HAAROPMIZATION

void VinisMatrixNormalization(double **mat, uint n, bool standard, bool invert = false);
void VinisNonStandardMatrixNormalization(double **matrix, uint n, bool invert = false);
void VinisNormalization(double *vec, uint n, bool invert = false);
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
void INT_Haar_Compression(interval *vec, int n, float percentage);
void INT_Haar_Matrix_Compression(interval **matrix, int n, float percentage);

#ifdef HAAROPMIZATION

void INT_VinisMatrixNormalization(interval **mat, uint n, bool standard, bool invert = false);
void INT_VinisNonStandardMatrixNormalization(interval **matrix, uint n, bool invert = false);
void INT_VinisNormalization(interval *vec, uint n, bool invert = false);
void INT_VinisStandardMatrixNormalization(interval **mat, uint n, bool invert = false);

#endif /* HAAROPMIZATION */

real INT_diameter(interval x);
real INT_error(interval *x, int n);
real INT_error(interval **x, int linhas, int colunas);

#endif /* INTHAAR */

#endif /* INTHAAR_H */