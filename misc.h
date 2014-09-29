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

void startTimeCounter();
double getTimeCounter();
void escrever_imagem(char *arquivo, double **matriz, ImageInfo imgInfo);
double** carregar_imagem(char *arquivo, ImageInfo *imageInfo);
void gnuplot_dat(const char *filename, double *x, double *y, int n);
void gnuplot_dat_Vdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_Wdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_VWdecomposition(const char *file1, const char *file2, double x1, double x2, double *v, int n, int levels, bool normal);
void data_analysis(double *data, uint n, DataAnalysis *analysis);

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