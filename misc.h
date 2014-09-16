#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <fstream>
#include "int-dwts.h"

struct ImageInfo
{
	int x, y;
	char magic[3];
};

void startTimeCounter();
double getTimeCounter();
void escrever_imagem(char *arquivo, double **matriz, struct ImageInfo imgInfo);
struct ImageInfo carregar_imagem(char *arquivo, double **data);
void gnuplot_dat(const char *filename, double *x, double *y, int n);
void gnuplot_dat_Vdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_Wdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal);
void gnuplot_dat_VWdecomposition(const char *file1, const char *file2, double x1, double x2, double *v, int n, int levels, bool normal);
