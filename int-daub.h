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
void Daub_DecompositionStep(T *v, uint n, bool normal, bool optimalFilters = true)
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

	T a, b, c, d, t1, t2, *offset;

	if (optimalFilters)
	{
		offset = v;

		for (i = 0, j = 0; j < n - 3; j += 2, i++)
		{
			a = *(offset++);
			b = *(offset++);
			c = *(offset++);
			d = *offset;

			t1 = a - c;
			t2 = b - d;

			result[i] = (a + d + v3*(b + c) + sqrt3*(t1 + t2)) / denom;
			result[i+half] = (a - d + v3*(c - b) + sqrt3*(t2 - t1)) / denom;

			offset--;
		}

		offset = v + n - 2;
		a = *(offset++);
		b = *offset;
		c = *v;
		d = *(v + 1);

		t1 = a - c;
		t2 = b - d;

		result[i] = (a + d + v3*(b + c) + sqrt3*(t1 + t2)) / denom;
		result[i+half] = (a - d + v3*(c - b) + sqrt3*(t2 - t1)) / denom;
	}
	else
	{
		denom = T(4) * sqrt(T(2));
		T h0 = (T(1) + sqrt3) / denom;
		T h1 = (T(3) + sqrt3) / denom;
		T h2 = (T(3) - sqrt3) / denom;
		T h3 = (T(1) - sqrt3) / denom;
		T g0 = h3;
		T g1 = -h2;
		T g2 = h1;
		T g3 = -h0;

		offset = v;

		for (i = 0, j = 0; j < n - 3; j += 2, i++)
		{
			a = *(offset++);
			b = *(offset++);
			c = *(offset++);
			d = *offset;

			t1 = a - c;
			t2 = b - d;

			result[i] = a*h0 + b*h1 + c*h2 + d*h3;
			result[i+half] = a*g0 + b*g1 + c*g2 + d*g3;

			offset--;
		}

		offset = v + n - 2;
		a = *(offset++);
		b = *offset;
		c = *v;
		d = *(v + 1);

		t1 = a - c;
		t2 = b - d;

		result[i] = a*h0 + b*h1 + c*h2 + d*h3;
		result[i+half] = a*g0 + b*g1 + c*g2 + d*g3;
	}

	for (i = 0; i < n; i++)
	{
		v[i] = result[i];
	}

	delete[] result;
}

template <typename T>
void Daub_Decomposition(T *vec, uint n, bool normal, bool optimalFilters = true)
{
	for (uint size = n; size >= 4; size >>= 1)
	{
		Daub_DecompositionStep(vec, size, normal, optimalFilters);
	}
}

template <typename T>
void Daub_CompositionStep(T *v, uint n, bool normal, bool optimalFilters = true)
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

	T a, b, c, d, t1, t2;
	T *offset1, *offset2;

	if (optimalFilters)
	{
		a = *(v + half - 1);
		b = *(v + n - 1);
		c = *v;
		d = *(v + half);

		t1 = c - a;
		t2 = b - d;

		result[0] = (c + d + v3*(a + b) + sqrt3*(t1 + t2)) / denom;
		result[1] = (a - b + v3*(c - d) + sqrt3*(t1 - t2)) / denom;

		for (i = 0, j = 2; i < half - 1; i++)
		{
			offset1 = v + i;
			offset2 = v + half + i;

			a = *(offset1++);
			b = *(offset2++);
			c = *(offset1++);
			d = *(offset2++);

			t1 = c - a;
			t2 = b - d;

			result[j++] = (c + d + v3*(a + b) + sqrt3*(t1 + t2)) / denom;
			result[j++] = (a - b + v3*(c - d) + sqrt3*(t1 - t2)) / denom;
		}
	}
	else
	{
		denom = T(4) * sqrt(T(2));
		T h0 = (T(1) + sqrt3) / denom;
		T h1 = (T(3) + sqrt3) / denom;
		T h2 = (T(3) - sqrt3) / denom;
		T h3 = (T(1) - sqrt3) / denom;
		T g0 = h3;
		T g1 = -h2;
		T g2 = h1;
		T g3 = -h0;

		T Ih0 = h2;
		T Ih1 = g2;
		T Ih2 = h0;
		T Ih3 = g0;
		T Ig0 = h3;
		T Ig1 = g3;
		T Ig2 = h1;
		T Ig3 = g1;

		a = *(v + half - 1);
		b = *(v + n - 1);
		c = *v;
		d = *(v + half);

		result[0] = a * Ih0 + b * Ih1 + c * Ih2 + d * Ih3;
		result[1] = a * Ig0 + b * Ig1 + c * Ig2 + d * Ig3;

		for (i = 0, j = 2; i < half - 1; i++)
		{
			offset1 = v + i;
			offset2 = v + half + i;

			a = *(offset1++);
			b = *(offset2++);
			c = *(offset1++);
			d = *(offset2++);

			result[j++] = a * Ih0 + b * Ih1 + c * Ih2 + d * Ih3;
			result[j++] = a * Ig0 + b * Ig1 + c * Ig2 + d * Ig3;
		}
	}

	for (i = 0; i < n; i++)
	{
		v[i] = result[i];
	}

	delete [] result;
}

template <typename T>
void Daub_Composition(T *vec, uint n, bool normal, bool optimalFilters = true)
{
	for (uint size = 4; size <= n; size <<= 1)
	{
		Daub_CompositionStep(vec, size, normal, optimalFilters);
	}
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
		
		int currentLevel = (int)(levels - level);

		factor = pow(T(2), T(currentLevel) / T(2));

		if (invert)
		for (uint i = start; i < end; i++)
			vec[i] *= factor;
		else
		for (uint i = start; i < end; i++)
			vec[i] /= factor;
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

		int currentLevel = (int)(levels - level);

		factor = pow(T(2), T(currentLevel) / T(2));

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
void Daub_StandardDecomposition(T **mat, uint rows, uint cols, bool normal, bool stepnorm = false, bool optimalFilters = true)
{
	T *temp_row = new T[cols];
	T *temp_col = new T[rows];

	for (uint i = 0; i < rows; i++)
	{
		for (uint j = 0; j < cols; j++)
			temp_row[j] = mat[i][j];

		Daub_Decomposition(temp_row, cols, normal, optimalFilters);

		for (uint j = 0; j < cols; j++)
			mat[i][j] = temp_row[j];
	}

	if (!normal && stepnorm && optimalFilters)
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

		Daub_Decomposition(temp_col, rows, normal, optimalFilters);

		for (uint j = 0; j < rows; j++)
			mat[j][i] = temp_col[j];
	}

	if (!normal && stepnorm && optimalFilters)
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
void Daub_StandardComposition(T **mat, uint rows, uint cols, bool normal, bool optimalFilters = true)
{
	T *temp_row = new T[cols];
	T *temp_col = new T[rows];

	for (uint i = 0; i < cols; i++)
	{
		for (uint j = 0; j < rows; j++)
			temp_col[j] = mat[j][i];

		Daub_Composition(temp_col, rows, normal, optimalFilters);

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

		Daub_Composition(temp_row, cols, normal, optimalFilters);

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
void Daub_NonStandardDecomposition(T **mat, uint rows, uint cols, bool normal, bool optimalFilters = true)
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

				Daub_DecompositionStep(temp_row, w, normal, optimalFilters);

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

				Daub_DecompositionStep(temp_col, h, normal, optimalFilters);

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
void Daub_NonStandardComposition(T **mat, uint rows, uint cols, bool normal, bool optimalFilters = true)
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

				Daub_CompositionStep(temp_col, r, normal, optimalFilters);

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

				Daub_CompositionStep(temp_row, c, normal, optimalFilters);

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
