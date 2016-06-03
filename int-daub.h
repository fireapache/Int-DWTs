#pragma once

#include <iostream>
#include <cmath>
#include "common.h"

using namespace std;

#ifndef WIN32
#include <interval.hpp>
using namespace cxsc;
#endif

template <typename T>
void Daub_Decomposition(T *vec, uint n, bool normal)
{
	for (uint size = n; size >= 4; size >>= 1)
	{
		Daub_DecompositionStep(vec, size, normal);
	}
}

template <typename T>
void Daub_DecompositionStep(T *v, uint n, bool normal)
{
	if (n < 4) return;

	T v3 = T(3);
	T sqrt3 = sqrt(T(3));
	T denom;
	uint i, j;
	const uint half = n >> 1;

	if (normal) denom = T(4) * sqrt(T(2));
	else denom = T(4);

	T* result = new T[n];

	for (i = 0, j = 0; j < n - 3; j += 2, i++)
	{
		result[i] = (v[j] + v[j+3] + v3*(v[j+1] + v[j+2]) + sqrt3*(v[j] + v[j+1] - v[j+2] - v[j+3])) / denom;
		result[i+half] = (v[j] - v[j+3] + v3*(v[j+2] - v[j+1]) + sqrt3*(v[j+1] - v[j] + v[j+2] - v[j+3])) / denom;
	}

	result[i] = (v[n-2] + v[1] + v3*(v[n-1] + v[0]) + sqrt3*(v[n-2] + v[n-1] - v[0] - v[1])) / denom;
	result[i+half] = (v[n-2] - v[1] + v3*(v[0] - v[n-1]) + sqrt3*(v[n-1] - v[n-2] + v[0] - v[1])) / denom;

	for (i = 0; i < n; i++)
	{
		v[i] = result[i];
	}

	delete[] result;
}

template <typename T>
void Daub_Composition(T *vec, uint n, bool normal)
{
	for (uint size = 4; size <= n; size <<= 1)
	{
		Daub_CompositionStep(vec, size, normal);
	}
}

template <typename T>
void Daub_CompositionStep(T *v, uint n, bool normal)
{
	if (n < 4) return;

	T v3 = T(3);
	T sqrt3 = sqrt(T(3));
	T denom;
	uint i, j;
	const uint half = n >> 1;
	const uint halfPls1 = half + 1;

	if (normal) denom = T(4) * sqrt(T(2));
	else denom = T(4);

	T* result = new T[n];

	result[0] = (v[0] + v[half] + v3*(v[half - 1] + v[n - 1]) + sqrt3*(v[n - 1] + v[0] - v[half - 1] - v[half])) / denom;
	result[1] = (v[half - 1] - v[n - 1] + v3*(v[0] - v[half]) + sqrt3*(v[0] + v[half] - v[half - 1] - v[n - 1])) / denom;

	for (i = 0, j = 2; i < half - 1; i++)
	{
		result[j++] = (v[i + 1] + v[i + halfPls1] + v3*(v[i] + v[i + half]) + sqrt3*(v[i + half] + v[i + 1] - v[i] - v[i + halfPls1])) / denom;
		result[j++] = (v[i] - v[i + half] + v3*(v[i + 1] - v[i + halfPls1]) + sqrt3*(v[i + 1] + v[i + halfPls1] - v[i] - v[i + half])) / denom;
	}

	for (i = 0; i < n; i++)
	{
		v[i] = result[i];
	}

	delete [] result;
}

template <typename T>
void Daub_Normalization(T *vec, uint n, bool invert = false)
{
	uint levels = (uint)log2(n);
	T factor;
	
	for (uint level = 1; level < levels; level++)
	{
		uint start;
		
		if (level == 1) start = 0;
		else start = (uint)pow(2.0, (double)(level));
		
		uint end = (uint)pow(2.0, (double)(level + 1));
		
		factor = pow(T(2), T(levels - level) / T(2));

		if (invert)
		for (uint i = start; i < end; i++)
			vec[i] *= factor;
		else
		for (uint i = start; i < end; i++)
			vec[i] /= factor;
	}
}

template <typename T>
void Daub_StandardDecomposition(T *mat, uint rows, uint cols, bool normal, bool stepnorm = false)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	for (uint i = 0; i < rows; i++)
	{
		for (uint j = 0; j < cols; j++)
			temp_row[j] = mat[i][j];

		Daub_Decomposition(temp_row, cols, normal);

		for (uint j = 0; j < cols; j++)
			mat[i][j] = temp_row[j];
	}

	if (!normal && stepnorm)
	{
		Daub_StandardStepNormalization(mat, rows, cols, true);
	}

	//if (normal) cout << "Normalized ";
	//cout << "Decomposition (rows):" << endl;
	//cout << endl;
	//printMatrix(mat, rows);
	//cout << endl;

	for (uint i = 0; i < cols; i++)
	{
		for (uint j = 0; j < rows; j++)
			temp_col[j] = mat[j][i];

		Daub_Decomposition(temp_col, rows, normal);

		for (uint j = 0; j < rows; j++)
			mat[j][i] = temp_col[j];
	}

	if (!normal && stepnorm)
	{
		Daub_StandardStepNormalization(mat, rows, cols, false);
	}

	//if (normal) cout << "Normalized ";
	//cout << "Decomposition (columns):" << endl;
	//cout << endl;
	//printMatrix(mat, rows);
	//cout << endl;

	delete[] temp_row;
	delete[] temp_col;
}

template <typename T>
void Daub_StandardComposition(T *mat, uint rows, uint cols, bool normal)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	for (uint i = 0; i < cols; i++)
	{
		for (uint j = 0; j < rows; j++)
			temp_col[j] = mat[j][i];

		Daub_Composition(temp_col, rows, normal);

		for (uint j = 0; j < rows; j++)
			mat[j][i] = temp_col[j];
	}

	//if (normal) cout << "Normalized ";
	//cout << "Composition (columns):" << endl;
	//cout << endl;
	//printMatrix(mat, rows);
	//cout << endl;

	for (uint i = 0; i < rows; i++)
	{
		for (uint j = 0; j < cols; j++)
			temp_row[j] = mat[i][j];

		Daub_Composition(temp_row, cols, normal);

		for (uint j = 0; j < cols; j++)
			mat[i][j] = temp_row[j];
	}

	//if (normal) cout << "Normalized ";
	//cout << "Composition (rows):" << endl;
	//cout << endl;
	//printMatrix(mat, rows);
	//cout << endl;

	delete[] temp_row;
	delete[] temp_col;
}

template <typename T>
void Daub_StandardNormalization(T **mat, uint n, bool invert = false)
{
	uint levels = (uint)log2(n);
	T factor;

	for (uint levelR = 1; levelR < levels; levelR++)
	{
		uint startR;

		if (levelR == 1)	startR = 0;
		else        		startR = (uint)pow(2.0, (double)(levelR));

		uint endR = (uint)pow(2.0, (double)(levelR + 1));

		for (uint row = startR; row < endR; row++)
		{
			for (uint levelC = 1; levelC < levels; levelC++)
			{
				uint startC;

				if (levelC == 1)	startC = 0;
				else        		startC = (uint)pow(2.0, (double)(levelC));

				uint endC = (uint)pow(2.0, (double)(levelC + 1));

				for (uint c = startC; c < endC; c++)
				{
					uint levelSum = (2 * levels - levelR - levelC);

					if (levelSum & 1)
					{
						if (levelSum > 2)	factor = T((int)(levelSum)-1) * sqrt(T(2));
						else				factor = sqrt(T(2));
					}
					else					factor = T((int)(levelSum));

					if (factor != T(0))
					{
						if (invert)	mat[row][c] *= factor;
						else 		mat[row][c] /= factor;
					}
				}
			}
		}
	}
}

template <typename T>
void Daub_StandardStepNormalization(T **mat, uint rows, uint cols, bool horizontal, bool invert = false)
{
	uint levels = (uint)log2(rows);
	T factor;

	for (uint level = 1; level < levels; level++)
	{
		uint start;

		if (level == 1) start = 0;
		else start = (uint)pow(2.0, (double)(level));

		uint end = (uint)pow(2.0, (double)(level + 1));

		factor = pow(T(2), T(levels - level) / T(2));

		if (invert)
		{
			if (horizontal)
			{
				for (uint i = 0; i < rows; i++)
				{
					for (uint j = start; j < end; j++)
					{
						mat[i][j] *= factor;
					}
				}
			}
			else
			{
				for (uint i = start; i < end; i++)
				{
					for (uint j = 0; j < cols; j++)
					{
						mat[i][j] *= factor;
					}
				}
			}
		}
		else
		{
			if (horizontal)
			{
				for (uint i = 0; i < rows; i++)
				{
					for (uint j = start; j < end; j++)
					{
						mat[i][j] /= factor;
					}
				}
			}
			else
			{
				for (uint i = start; i < end; i++)
				{
					for (uint j = 0; j < cols; j++)
					{
						mat[i][j] /= factor;
					}
				}
			}
		}
	}
}

template <typename T>
void Daub_NonStandardDecomposition(T **mat, uint rows, uint cols, bool normal)
{
	uint h = rows, w = cols;
	T *temp_row = new T[cols];
	T *temp_col = new T[rows];

	while (w >= 4 || h >= 4)
	{
		if (w >= 4)
		{
			for (uint i = 0; i < h; i++)
			{
				for (uint j = 0; j < w; j++)
					temp_row[j] = mat[i][j];

				Daub_DecompositionStep(temp_row, w, normal);

				for (uint j = 0; j < w; j++)
					mat[i][j] = temp_row[j];
			}
		}

		if (h >= 4)
		{
			for (uint i = 0; i < w; i++)
			{
				for (uint j = 0; j < h; j++)
					temp_col[j] = mat[j][i];

				Daub_DecompositionStep(temp_col, h, normal);

				for (uint j = 0; j < h; j++)
					mat[j][i] = temp_col[j];
			}
		}

		if (w >= 4) w /= 2;
		if (h >= 4) h /= 2;
	}

	delete [] temp_row;
	delete [] temp_col;
}

template <typename T>
void Daub_NonStandardComposition(T **mat, uint rows, uint cols, bool normal)
{
	uint r = 4, c = 4;
	T *temp_row = new T[cols];
	T *temp_col = new T[rows];

	while (c <= cols || r <= rows)
	{
		if (r <= rows)
		{
			for (uint i = 0; i < c; i++)
			{
				for (uint j = 0; j < rows; j++)
					temp_col[j] = mat[j][i];

				Daub_CompositionStep(temp_col, r, normal);

				for (uint j = 0; j < rows; j++)
					mat[j][i] = temp_col[j];
			}
		}

		if (c <= cols)
		{
			for (uint i = 0; i < r; i++)
			{
				for (uint j = 0; j < cols; j++)
					temp_row[j] = mat[i][j];

				Daub_CompositionStep(temp_row, c, normal);

				for (uint j = 0; j < cols; j++)
					mat[i][j] = temp_row[j];
			}
		}

		if (c <= cols) c *= 2;
		if (r <= rows) r *= 2;
	}

	delete [] temp_row;
	delete [] temp_col;
}

template <typename T>
void Daub_NonStandardNormalization(T **mat, uint n, bool invert = false)
{
	T factor;
	uint start, end;

	uint levels = (uint)log2((double)n);

	for (uint i = 1; i < levels; i++)
	{
		if (i == 1) start = 0;
		else start = (uint)pow(2.0, (double)i);
		end = (uint)pow(2.0, (double)i + 1);

		factor = pow(T(2.0), T(levels - i));

		if (i == 1)
		{
			if (invert)
			{
				for (uint l = 0; l < 4; l++)
				for (uint c = 0; c < 4; c++)
					mat[l][c] *= factor;
			}
			else
			{
				for (uint l = 0; l < 4; l++)
				for (uint c = 0; c < 4; c++)
					mat[l][c] /= factor;
			}

			continue;
		}

		if (invert)
		{
			for (uint l = 0; l < end / 2; l++)
			for (uint c = start; c < end; c++)
				mat[l][c] *= factor;
		}
		else
		{
			for (uint l = 0; l < end / 2; l++)
			for (uint c = start; c < end; c++)
				mat[l][c] /= factor;
		}

		if (invert)
		{
			for (uint l = start; l < end; l++)
			for (uint c = 0; c < end; c++)
				mat[l][c] *= factor;
		}
		else
		{
			for (uint l = start; l < end; l++)
			for (uint c = 0; c < end; c++)
				mat[l][c] /= factor;
		}

	}
}
