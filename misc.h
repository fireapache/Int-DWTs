#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

struct ImageInfo
{
	int x, y;
	char magic[3];
};

void startTimeCounter();
double getTimeCounter();
void escrever_imagem(char *arquivo, double **matriz);
struct ImageInfo carregar_imagem(char *arquivo, double **data);