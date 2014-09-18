#include "misc.h"

#ifdef WIN32
double pcFreq = 0.0;
__int64 counterStart = 0;
#else
timeval tCounter, tTime;
#endif

double** carregar_imagem(char *arquivo, ImageInfo *imageInfo)
{
    char auxc[5];
    int l, c, aux, temp[3];
    double **data;
    FILE *imagem;
    ImageInfo imgInfo;

    if ((imagem = fopen(arquivo, "r")) == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", arquivo);
        return NULL;
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

    data = new double*[imgInfo.x];

    for (int i = 0; i < imgInfo.x; i++)
    {
    	data[i] = new double[imgInfo.x];
    }
    
    for (l = 0; l < imgInfo.y; l++)
    for (c = 0; c < imgInfo.x; c++)
    {
        fscanf(imagem, "%d", &temp[0]);
        fscanf(imagem, "%d", &temp[1]);
        fscanf(imagem, "%d", &temp[2]);
        data[l][c] = ((double)temp[2] + (double)temp[1] + (double)temp[0]) / 3.0;
    }

    imageInfo->x = imgInfo.x;
    imageInfo->y = imgInfo.y;
    imageInfo->magic[0] = imgInfo.magic[0];
    imageInfo->magic[1] = imgInfo.magic[1];
    imageInfo->magic[2] = imgInfo.magic[2];

    return data;
}

void startTimeCounter()
{
#ifdef WIN32
	LARGE_INTEGER li;
	QueryPerformanceFrequency(&li);
	pcFreq = double(li.QuadPart);
	QueryPerformanceCounter(&li);
	counterStart = li.QuadPart;
#else
	gettimeofday(&tCounter, NULL);
#endif
}

double getTimeCounter()
{
	double returnTime;

#ifdef WIN32
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	returnTime = double(li.QuadPart - counterStart) / pcFreq;
#else
	gettimeofday(&tTime, NULL);
    double usecDiff = tTime.tv_usec - tCounter.tv_usec;
    if (usecDiff < 0) returnTime = (double)(usecDiff + 1000000) / 1000000.0;
    else returnTime = usecDiff / 1000000.0;
#endif

    return returnTime;
}

void escrever_imagem(char *arquivo, double **matriz, ImageInfo imgInfo)
{
	int data;
    ofstream out;

    out.open(arquivo, ios_base::trunc);

    out << imgInfo.magic << endl;
    out << imgInfo.x << " " << imgInfo.y << endl;
    out << "255" << endl;

    for (int l = 0; l < imgInfo.y; l++)
    {
        for (int c = 0; c < imgInfo.x; c++)
        {
			data = static_cast<int>(matriz[l][c]);
            out << data << " " << data << " " << data << endl;
        }
    }

    out.close();
}

void gnuplot_dat(const char *filename, double *x, double *y, int n)
{
    ofstream out;

    out.open(filename, ios_base::trunc);

    for (int i = 0; i < n; i++)
        out << x[i] << '\t' << y[i] << '\n';

    out.close();
}

void gnuplot_dat_Vdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal)
{
    ofstream out;
    double x, alpha;

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] / sqrt(float(n));

    for (int i = 0; i < levels; i++)
    {
        Haar_DecompositionStep(v, n, normal);
        n /= 2;
    }

    out.open(file, ios_base::trunc);

    for (int i = 0; i < n; i++)
    {
        alpha = (double)i / (double)(n);
        x = x1 * (1.0 - alpha) + x2 * alpha;

        out << x << '\t' << v[i] << '\n';
    }

    out.close();

    for (int i = levels; i > 0; i--)
    {
        Haar_CompositionStep(v, n, normal);
        n *= 2;
    }

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] * sqrt(float(n));
}

void gnuplot_dat_Wdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal)
{
    ofstream out;
    double x, alpha;

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] / sqrt(float(n));

    for (int i = 0; i < levels; i++)
    {
        Haar_DecompositionStep(v, n, normal);
        n /= 2;
    }

    out.open(file, ios_base::trunc);

    for (int i = 0; i < n; i++)
    {
        alpha = (double)i / (double)(n);
        x = x1 * (1.0 - alpha) + x2 * alpha;
        x += (x2 - x1) / (2 * n);
        out << x << '\t' << v[i + n] << '\n';
    }

    out.close();

    for (int i = levels; i > 0; i--)
    {
        Haar_CompositionStep(v, n, normal);
        n *= 2;
    }

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] * sqrt(float(n));
}

void gnuplot_dat_VWdecomposition(const char *file1, const char *file2, double x1, double x2, double *v, int n, int levels, bool normal)
{
    ofstream out1, out2;
    double x, alpha;

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] / sqrt(float(n));

    for (int i = 0; i < levels; i++)
    {
        Haar_DecompositionStep(v, n, normal);
        n /= 2;
    }

    out1.open(file1, ios_base::trunc);
    out2.open(file2, ios_base::trunc);

    for (int i = 0; i < n; i++)
    {
        alpha = (double)i / (double)(n);
        x = x1 * (1.0 - alpha) + x2 * alpha;

        out1 << x << '\t' << v[i] << '\n';
        out2 << x << '\t' << v[i + n] << '\n';
    }

    out1.close();
    out2.close();

    for (int i = levels; i > 0; i--)
    {
        Haar_CompositionStep(v, n, normal);
        n *= 2;
    }

    if (normal)
    for (int i = 0; i < n; i++)
        v[i] = v[i] * sqrt(float(n));
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