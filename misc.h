#ifndef MISC_H
#define MISC_H

#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <fstream>
#include <iomanip>
#include "int-dwts.h"

typedef struct ImageInfo
{
	int x, y;
	char magic[3];
} ImageInfo;

typedef struct DataAnalysis
{
	double mean;
	double deviation;
} DataAnalysis;

template <typename T>
struct EucMSE_data
{
	T mse;
	T euc;
};

template <typename T>
struct ImageQuality
{
	T mse;
	T euc;
	T psnr;
};

typedef struct TimeMesurement
{
	double mean;
	double dev;
	double stdDev;
} TimeMesurement;

template <typename T>
struct MatrixMap
{
	uint x;
	uint y;
	T value;
};

template <typename T>
void matrixMapping(T **matrix, MatrixMap<T> *map, uint n)
{
	for (uint i = 0; i < n; i++)
	{
		for (uint j = 0; j < n; j++)
		{
			map[(i*n) + j].x = i;
			map[(i*n) + j].y = j;
			map[(i*n) + j].value = matrix[i][j];
		}
	}
}

template <typename T>
void matrixRemapping(T **matrix, MatrixMap<T> *map, uint n)
{
	for (uint i = 0; i < n*n; i++)
	{
		matrix[map[i].x][map[i].y] = map[i].value;
	}
}

template <typename T>
void merging(MatrixMap<T> *a, MatrixMap<T> *b, uint low, uint mid, uint high)
{
	uint l1, l2, i;

	for (l1 = low, l2 = mid + 1, i = low; l1 <= mid && l2 <= high; i++)
	{
		if (a[l1].value <= a[l2].value)
		{
			b[i] = a[l1++];
		}
		else
		{
			b[i] = a[l2++];
		}
	}

	while (l1 <= mid)
	{
		b[i++] = a[l1++];
	}

	while (l2 <= high)
	{
		b[i++] = a[l2++];
	}

	for (i = low; i <= high; i++)
	{
		a[i] = b[i];
	}
}

template <typename T>
void sort(MatrixMap<T> *a, MatrixMap<T> *b, uint low, uint high)
{
	uint mid;

	if (low < high)
	{
		mid = (low + high) / 2;
		sort(a, b, low, mid);
		sort(a, b, mid + 1, high);
		merging(a, b, low, mid, high);
	}
	else
	{
		return;
	}
}

template <typename T>
T getMinValue(MatrixMap<T> *v, uint n)
{
	if (!v) return T(0);

	T minValue = v[0].value;

	for (uint i = 0; i < n; i++)
	{
		if (v[i].value < minValue)
		{
			minValue = v[i].value;
		}
	}

	return minValue;
}

template <typename T>
T getMaxValue(MatrixMap<T> *v, uint n)
{
	if (!v) return T(0);

	T maxValue = v[0].value;

	for (uint i = 0; i < n; i++)
	{
		if (v[i].value > maxValue)
		{
			maxValue = v[i].value;
		}
	}

	return maxValue;
}

template <typename T>
T getThreshold(MatrixMap<T> *v, uint n, double ratio)
{
	T threshold;
	T error;
	T totalEnergy = T(0);
	T tMin = getMinValue(v, n);
	T tMax = getMaxValue(v, n);

	if (ratio > 1.0) ratio = 1.0;
	else if (ratio < 0.0) ratio = 0.0;

	for (uint i = 0; i < n; i++)
	{
		totalEnergy += pow(static_cast<double>(v[i].value), 2.0);
	}

	error = totalEnergy * ratio;

	T temp;

	do
	{
		threshold = (tMax + tMin) / 2;
		totalEnergy = T(0);

		for (uint i = 0; i < n; i++)
		{
			if (v[i].value < threshold)
			{
				totalEnergy += pow(static_cast<double>(v[i].value), 2.0);
			}
			else break;
		}

		if (totalEnergy < error)
		{
			tMin = threshold;
		}
		else
		{
			tMax = threshold;
		}

		temp = tMax - tMin;

	} while ((tMax - tMin) > 0.001);

	return threshold;
}

template <typename T>
void compress(MatrixMap<T> *v, uint n, double ratio)
{
	T threshold = getThreshold(v, n, ratio);

	for (uint i = 0; i < n; i++)
	{
		if (v[i].value < threshold)
		{
			v[i].value = T(0);
		}
		else break;
	}

}

void startTimeCounter();
double getTimeCounter();
void escrever_imagem(char *arquivo, double **matriz, ImageInfo imgInfo);
double** carregar_imagem(char *arquivo, ImageInfo *imageInfo);
void gnuplot_dat(const char *filename, double *x, double *y, int n);
void gnuplot_dat_Vdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_Wdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_VWdecomposition(const char *file1, const char *file2, double x1, double x2, double *v, int n, int levels, bool normal);
void data_analysis(double *data, uint n, DataAnalysis *analysis);
void data_analysis(double **data, uint n, DataAnalysis *analysis);

TimeMesurement runTimeMesurement(double *times, uint n);

bool isPowerOfTwo(unsigned int x);

template <typename T>
T* matrixToVector(T **mat, uint n)
{
	T *vec = new T[n*n];

	for (uint i = 0; i < n; i++)
	{
		for (uint j = 0; j < n; j++)
		{
			vec[i*n + j] = mat[i][j];
		}
	}

	return vec;
}

template <typename T>
T** genRandomMatrix(uint rows, uint cols, uint max)
{
	T **result;
	double r;

	result = new T*[rows];

	for (uint i = 0; i < rows; ++i)
	{
		result[i] = new T[cols];

		for (uint j = 0; j < cols; ++j)
		{
			r = rand() % max;
			result[i][j] = T(r);
		}

	}

	return result;

}

template <typename T>
T*** genRandomMatrices(uint mats, uint rows, uint cols, uint max)
{
	T ***result;
	double r;

	result = new T**[mats];

	for (uint i = 0; i < mats; ++i)
	{
		result[i] = new T*[rows];

		for (uint j = 0; j < rows; ++j)
		{
			result[i][j] = new T[cols];
	
			for (uint w = 0; w < cols; ++w)
			{
				r = rand() % max;
				result[i][j][w] = T(r);
			}
		}
	}

	return result;

}

template <typename T>
void printVectors(T **vecs, uint n1, uint n2, uint space, bool row = true)
{
	if (n1 <= 0 || n2 <= 0 || space <= 0) return;

	if (row)
	for (uint i = 0; i < n2; ++i)
	{
		cout << vecs[i][0];

		for (uint j = 1; j < n1; ++j)
		{
			cout << setw(space) << vecs[i][j];
		}

		cout << endl;
	}
	else
	for (uint i = 0; i < n1; ++i)
	{
		cout << vecs[0][i];

		for (uint j = 1; j < n1; ++j)
		{
			cout << setw(space) << vecs[j][i];
		}

		cout << endl;
	}
}

template <typename T>
void printMatrices(T ***mats,  uint n1, uint n2, uint n3, uint space)
{
	if (n1 <= 0 || n2 <= 0 || n3 <= 0 || space <= 0) return;

	for (uint i = 0; i < n3; ++i)
	{
		printVectors(mats[i], n1, n2, space);
		cout << "====================" << endl;
	}
}

template <typename T>
void copyVector(T *v1, T *v2, uint n)
{
	for (uint i = 0; i < n; i++)
	{
		v2[i] = v1[i];
	}
}

template <typename T>
void copyMatrix(T **m1, T **m2, uint n)
{
	for (uint i = 0; i < n; i++)
	for (uint j = 0; j < n; j++)
	{
		m2[i][j] = m1[i][j];
	}
}

template <typename T>
void deleteMatrix(T **matrix, uint n)
{
	for (uint i = 0; i < n; i++)
	{
		delete [] matrix[i];
	}

	delete [] matrix;
}

template <typename T>
void deleteMatrices(T ***matrices, uint nMatrices, uint n)
{
	for (uint i = 0; i < nMatrices; ++i)
	{
		deleteMatrix<T>(matrices[i], n);
	}

	delete [] matrices;
}

template <typename T>
T MSE(T **m1, T **m2, uint n)
{
	T result = T();

	for (uint i = 0; i < n; ++i)
	for (uint j = 0; j < n; ++j)
	{
		result += pow(m1[i][j] - m2[i][j], T(2.0));
	}

	return result / pow(n, T(2.0));
}

template <typename T>
EucMSE_data<T> EucMSE(T *vec1, T *vec2, uint n)
{
	EucMSE_data<T> result;

	result.mse = T(0.0);
	result.euc = T(0.0);

	for (uint i = 0; i < n; ++i)
	{
		result.mse += pow(vec1[i] - vec2[i], T(2.0));
		result.euc += pow(vec1[i] - vec2[i], T(2.0));
	}

	result.mse = result.mse / n;
	result.euc = sqrt(result.euc);

	return result;
}

template <typename T>
EucMSE_data<T> EucMSE(T **m1, T **m2, uint n)
{
	EucMSE_data<T> result;	

	result.mse = T(0.0);
	result.euc = T(0.0);

	for (uint i = 0; i < n; ++i)
	for (uint j = 0; j < n; ++j)
	{
		result.mse += pow(m1[i][j] - m2[i][j], T(2.0));
		result.euc += pow(m1[i][j] - m2[i][j], T(2.0));
	}

	result.mse = result.mse / pow(n, T(2.0));
	result.euc = sqrt(result.euc);

	return result;
}

EucMSE_data<interval> INT_EucMSE(interval **m1, interval **m2,  int n);

template <typename T>
T PSNR(T mse, T max)
{
	T result = T(0.0);

	if (mse != T(0.0))
	{
		result = 10.0 * log10((max * max) / mse);
	}
	else result = T(0.0);

	return result;
}

interval INT_PSNR(interval mse, interval max);

template <typename T>
ImageQuality<T> imageQuality(T *vec1, T *vec2, int max, uint n)
{
	ImageQuality<T> result;
    EucMSE_data<T> data;

    data = EucMSE<T>(vec1, vec2, n);

    result.mse = data.mse;
    result.euc = data.euc;
    result.psnr = PSNR<T>(data.mse, T(max));

    return result;
}

template <typename T>
T minValue(T *vec, uint n)
{
	T ref = vec[0];

	for (uint i = 0; i < n; i++)
	{
		if (ref > vec[i])
		{
			ref = vec[i];
		}
	}

	return ref;
}

template <typename T>
T maxValue(T *vec, uint n)
{
	T ref = vec[0];

	for (uint i = 0; i < n; i++)
	{
		if (ref < vec[i])
		{
			ref = vec[i];
		}
	}

	return ref;
}

template <typename T>
T maxValue(T **mat, uint n, uint m)
{
	T ref = mat[0][0];

	for (uint i = 0; i < n; i++)
	for (uint j = 0; j < m; j++)
	{
		if (ref < mat[i][j])
		{
			ref = mat[i][j];
		}
	}

	return ref;
}

template <typename T>
ImageQuality<T> imageQuality(T **m1, T **m2, int max, uint n)
{
	ImageQuality<T> result;
    EucMSE_data<T> data;

    data = EucMSE<T>(m1, m2, n);

    result.mse = data.mse;
    result.euc = data.euc;
    result.psnr = PSNR<T>(data.mse, T(max));

    return result;
}

template <typename T>
void printVector(T *vec, int n)
{
	int i;

	for (i = 0; i < n; i++)
	{
		cout << vec[i] << endl;
	}

	cout << endl;
}

template <typename T>
void printMatrix(T **mat, int n)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			cout << mat[i][j] << '\t';
		}

		cout << endl;
	}

	cout << endl;
}

// Speedup of A over B.
template <typename T>
T speedupCalc(T A, T B)
{
	return T(100) * (T(1) - A / B);
}

// Relative gain of A over B
template <typename T>
T relativeGain(T A, T B)
{
	return (B - A) / B * T(100);
}

#endif /* MISC_H */