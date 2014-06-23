#include "misc.h"

struct ImageInfo carregar_imagem(char *arquivo, double **data)
{
    char auxc[5];
    int l, c, aux, temp[3];
    FILE *imagem;
    struct ImageInfo imgInfo;

    if ((imagem = fopen(arquivo, "r")) == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", arquivo);
        return imgInfo;
    }

    fscanf(imagem, "%s", imgInfo.magic);
    fgetc(imagem);

    fscanf(imagem,"%s", auxc);
    if (auxc[0]=='#')
    {
        while(fgetc(imagem) != '\n');
        fscanf (imagem,"%d", &imgInfo.x);
    }
    else imgInfo.x = atoi(auxc);

	fscanf(imagem,"%d", &imgInfo.y);
    fscanf(imagem, "%d", &aux);

    if ((data = (double **)malloc(sizeof(double *)*imgInfo.y)) == NULL)
    {
        printf("Erro na alocação de memória!1\n");
        return imgInfo;
    }
    for (aux = 0; aux < imgInfo.x; aux++)
    if ((data[aux] = (double *)malloc(sizeof(double)*imgInfo.x)) == NULL)
    {
        printf("Erro na alocação de memória!\n");
        return imgInfo;
    }
    
    for (l = 0; l < imgInfo.y; l++)
    for (c = 0; c < imgInfo.x; c++)
    {
        fscanf(imagem, "%d", &temp[0]);
        fscanf(imagem, "%d", &temp[1]);
        fscanf(imagem, "%d", &temp[2]);
        data[l][c] = ((double)temp[2] + (double)temp[1] + (double)temp[0]) / 3.0;
    }

    return imgInfo;
}

void escrever_imagem(char *arquivo, double **matriz, struct ImageInfo imgInfo)
{
    FILE *imagem;
    int l, c, auxi;

    if ((imagem = fopen(arquivo, "w")) == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", arquivo);
        return;
    }

    fwrite(&imgInfo.magic[0], 1, 1, imagem);
    fwrite(&imgInfo.magic[1], 1, 1, imagem);
    fprintf(imagem, "\n%d %d\n255\n", imgInfo.x, imgInfo.y);

    for (l = 0; l < imgInfo.y; l++)
    {
        for (c = 0; c < imgInfo.x; c++)
        {
			auxi = static_cast<int>(matriz[l][c]);
            fprintf(imagem, "%d %d %d ", auxi, auxi, auxi);
        }

        fprintf(imagem, "\n");
    }

    fclose(imagem);
}

/*void escrever_imagem_from_greyscale(char *arquivo, double **matriz)
{
    FILE *imagem;
    //unsigned char auxc;
    int l, c, auxi, auxc;

    if ((imagem = fopen(arquivo, "w")) == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", arquivo);
        return;
    }

    fwrite(&magic[0], 1, 1, imagem);
    fwrite(&magic[1], 1, 1, imagem);
    fprintf(imagem, "\n%d %d\n255\n", x, y);
    for (l = 0; l < y; l++)
    {
        for (c = 0; c < x; c++)
        {
			auxi = (int)matriz[l][c];
			
            fprintf(imagem, "%d\n", auxi & 0xFF);
            fprintf(imagem, "%d\n", auxi & 0xFF);
            fprintf(imagem, "%d\n", auxi & 0xFF);
        }
    }

    fclose(imagem);
}

// Transforma as cores da imagem em tons de cinza.
double **to_greyscale(double **matrix, int x, int y)
{
	int aux[3], aux2, aux3;
	double **retorno = new double*[y];

	for (int m = 0; m < y; m++)
	retorno[m] = new double[x];
	
	for (int i = 0; i < x; i++)
	for (int j = 0; j < y; j++)
	{
		aux2 = (int)matrix[i][j];

		aux[0] = aux2 & 0xFF;
		aux2 /= 0x100;
		aux[1] = aux2 & 0xFF;
		aux2 /= 0x100;
		aux[2] = aux2 & 0xFF;

		aux2 = (aux[0] + aux[1] + aux[2]) / 3;

		//aux3 = aux2 * 0x10000 + aux2 * 0x100 + aux2;

		//matrix[i][j] = (double)aux3;
		retorno[i][j] = (double)aux2;
	}
	
	return retorno;
}

double **from_greyscale(double **matrix, int x, int y)
{
	int aux1, aux2;
	double **retorno = new double*[y];

	for (int m = 0; m < y; m++)
	retorno[m] = new double[x];
	
	for (int i = 0; i < x; i++)
	for (int j = 0; j < y; j++)
	{
		aux1 = (int)matrix[i][j];
		
		aux2 = aux1 * 0x100;
		aux2 = aux2 + (aux1 & 0xFF);
		aux2 = aux2 * 0x100;
		aux2 = aux2 + (aux1 & 0xFF);
		
		retorno[i][j] = (double)aux2;
	}
	
	return retorno;
}*/