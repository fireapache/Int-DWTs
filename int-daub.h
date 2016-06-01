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
