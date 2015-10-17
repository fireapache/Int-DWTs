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

template <typename T>
void printVectors(T **vecs, uint n1, uint n2, uint space, bool row = true)
{
	if (n1 <= 0 || n2 <= 0 || space <= 0) return;

	if (row)
	for (int i = 0; i < n2; ++i)
	{
		cout << vecs[i][0];

		for (int j = 1; j < n1; ++j)
		{
			cout << setw(space) << vecs[i][j];
		}

		cout << endl;
	}
	else
	for (int i = 0; i < n1; ++i)
	{
		cout << vecs[0][i];

		for (int j = 1; j < n1; ++j)
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

	for (int i = 0; i < n3; ++i)
	{
		printVectors(mats[i], n1, n2, space);
		cout << "====================" << endl;
	}
}

template <typename T>
void copyVector(T *v1, T *v2, uint n)
{
	for (int i = 0; i < n; i++)
	{
		v2[i] = v1[i];
	}
}

template <typename T>
void copyMatrix(T **m1, T **m2, uint n)
{
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	{
		m2[i][j] = m1[i][j];
	}
}

template <typename T>
void deleteMatrix(T **matrix, uint n)
{
	for (int i = 0; i < n; i++)
	{
		delete [] matrix[i];
	}

	delete [] matrix;
}

template <typename T>
void deleteMatrices(T ***matrices, uint nMatrices, uint n)
{
	for (int i = 0; i < nMatrices; ++i)
	{
		deleteMatrix<T>(matrices[i], n);
	}
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

#endif /* MISC_H */