#include "int-haar.h"

void VinisNormalization(double *vec, UINT n)
{
	UINT levels = log2(n);
	
	for (UINT level = 0; level < levels; level++)
	{
		UINT inicio;
		
		if (level == 0) inicio = 0;
		else        	inicio = (UINT)pow(2.0, (double)(level));
		
		UINT fim    = (UINT)pow(2.0, (double)(level + 1));
		
		for (UINT i = inicio; i < fim; i++)
			vec[i] *= pow(2.0, - (double)level / 2.0);
	}
}

void INT_VinisNormalization(interval *vec, UINT n)
{
	UINT levels = (UINT)log2((double)(n));
	interval div;
	
	for (UINT level = 0; level < levels; level++)
	{
		UINT inicio;
		
		if (level == 0) inicio = 0;
		else			inicio = (UINT)pow(2.0, (double)(level));
		
		UINT fim    = (UINT)pow(2.0, (double)(level + 1));
		
		for (UINT i = inicio; i < fim; i++)
		{
			if (level & 1)
			{
				if (level > 2)	div = interval((int)(level) - 1) * sqrt(interval(2));
				else			div = sqrt(interval(2));
			}
			else				div = interval((int)(level));

			if (div != interval(0))
				vec[i] /= div;
		}
	}
}

void VinisStandardMatrixNormalization(double **mat, UINT n, bool invert)
{
	UINT levels = log2(n);
	
	for (UINT levelL = 0; levelL < levels; levelL++)
	{
		UINT inicioL;
		
		if (levelL == 0)	inicioL = 0;
		else        		inicioL = (UINT)pow(2.0, (double)(levelL));
		
		UINT fimL = (UINT)pow(2.0, (double)(levelL + 1));

		for (UINT l = inicioL; l < fimL; l++)
		{
			for (UINT levelC = 0; levelC < levels; levelC++)
			{
				UINT inicioC;
				
				if (levelC == 0)	inicioC = 0;
				else        		inicioC = (UINT)pow(2.0, (double)(levelC));
				
				UINT fimC = (UINT)pow(2.0, (double)(levelC + 1));

				if (invert)
				{
					for (UINT c = inicioC; c < fimC; c++)
						mat[l][c] *= pow(2.0, (double)(levelL + levelC) / 2.0);
				}
				else
				{
					for (UINT c = inicioC; c < fimC; c++)
						mat[l][c] /= pow(2.0, (double)(levelL + levelC) / 2.0);
				}
			}
		}
	}
}

void INT_VinisStandardMatrixNormalization(interval **mat, UINT n, bool invert)
{
	UINT levels = log2(n);
	interval div;
	
	for (UINT levelL = 0; levelL < levels; levelL++)
	{
		UINT inicioL;
		
		if (levelL == 0)	inicioL = 0;
		else        		inicioL = (UINT)pow(2.0, (double)(levelL));
		
		UINT fimL = (UINT)pow(2.0, (double)(levelL + 1));

		for (UINT l = inicioL; l < fimL; l++)
		{
			for (UINT levelC = 0; levelC < levels; levelC++)
			{
				UINT inicioC;
				
				if (levelC == 0)	inicioC = 0;
				else        		inicioC = (UINT)pow(2.0, (double)(levelC));
				
				UINT fimC = (UINT)pow(2.0, (double)(levelC + 1));

				for (UINT c = inicioC; c < fimC; c++)
				{
					UINT levelSum = (levelL + levelC);

					if (levelSum & 1)
					{
						if (levelSum > 2)	div = interval((int)(levelSum) - 1) * sqrt(interval(2));
						else				div = sqrt(interval(2));
					}
					else					div = interval((int)(levelSum));

					if (div != interval(0))
					{
						if (invert)	mat[l][c] *= div;
						else 		mat[l][c] /= div;
					}
				}
					//mat[l][c] *= pow(2.0, - (double)(levelL + levelC) / 2.0);
			}
		}
	}
}

void VinisNonStandardMatrixNormalization(double **matrix, UINT n, bool invert)
{
    double div;
    UINT start, limit;

    UINT level = (UINT)(log((double)n) / log(2.0)) - 1;

    if (level <= 0) return;

    for (UINT i = level; i > 0; i--)
    {
        div   = pow(2.0, (double)i);
        start = (UINT)div;
        limit = (UINT)pow(2.0, (double)i + 1);

        for (UINT l = 0; l < limit / 2; l++)
        {
            for (UINT c = start; c < limit; c++)
            {
                matrix[l][c] /= div;
            }
        }

        if (invert)
        {
        	for (UINT l = start; l < limit; l++)
        	for (UINT c = 0; c < limit; c++)
        		matrix[l][c] *= div;
        }
        else
        {
        	for (UINT l = start; l < limit; l++)
        	for (UINT c = 0; c < limit; c++)
        		matrix[l][c] /= div;
        }
        
    }
}

void INT_VinisNonStandardMatrixNormalization(interval **matrix, UINT n, bool invert)
{
	interval div;
    UINT start, limit;

    UINT level = (UINT)(log((double)n) / log(2.0)) - 1;

    if (level <= 0) return;

    for (UINT i = level; i > 0; i--)
    {
        div   = interval(pow(2.0, (double)i));
        start = (UINT)pow(2.0, (double)i);
        limit = (UINT)pow(2.0, (double)i + 1);

        for (UINT l = 0; l < limit / 2; l++)
        {
            for (UINT c = start; c < limit; c++)
            {
                matrix[l][c] /= div;
            }
        }

        if (invert)
        {
        	for (UINT l = start; l < limit; l++)
        	for (UINT c = 0; c < limit; c++)
        		matrix[l][c] *= div;
        }
        else
        {
        	for (UINT l = start; l < limit; l++)
        	for (UINT c = 0; c < limit; c++)
        		matrix[l][c] /= div;
        }
    }
}

void VinisMatrixNormalization(double **mat, UINT n, bool standard, bool invert)
{
	if (standard)	VinisStandardMatrixNormalization(mat, n, invert);
	else			VinisNonStandardMatrixNormalization(mat, n, invert);
}

void INT_VinisMatrixNormalization(interval **mat, UINT n, bool standard, bool invert)
{
	if (standard)	INT_VinisStandardMatrixNormalization(mat, n, invert);
	else			INT_VinisNonStandardMatrixNormalization(mat, n, invert);
}

// Passo de decomposição de tipo Haar de um vetor.
void Haar_DecompositionStep(double *vec, int n, bool normal)
{
	double div;
	double *vecp = new double[n];
	
	for (int i = 0; i < n; i++)
		vecp[i] = 0;
	
	if (normal) div = sqrt(2.0);
	else        div = 2.0;
	
	for (int i = 0; i < n/2; i++)
	{
		vecp[i]       = (vec[2*i] + vec[2*i + 1]) / div;
		vecp[i + n/2] = (vec[2*i] - vec[2*i + 1]) / div;
	}
	
	for (int i = 0; i < n; i++)
		vec[i] = vecp[i];
	
	delete [] vecp;
}

void INT_Haar_DecompositionStep(interval *vec, int n, bool normal)
{
	interval div;
	interval *vecp = new interval[n];
	
	for (int i = 0; i < n; i++)
		vecp[i] = interval(0);
	
	if (normal) div = sqrt(interval(2));
	else        div = interval(2);
	
	for (int i = 0; i < n/2; i++)
	{
		vecp[i]       = (vec[2*i] + vec[2*i + 1]) / div;
		vecp[i + n/2] = (vec[2*i] - vec[2*i + 1]) / div;
	}
	
	for (int i = 0; i < n; i++)
		vec[i] = vecp[i];
	
	delete [] vecp;
}

// Decomposição completa de tipo Haar de um vetor.
double **Haar_Decomposition_For_Graphs(double *vec, int n, bool normal)
{
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] / sqrt(float(n));

	double **_return = (double **)malloc(sizeof(double *)*2);

	_return[0] = (double *)malloc(sizeof(double)*(n - 1));
	_return[1] = (double *)malloc(sizeof(double)*(n - 1));
	
	while (n > 1)
	{
		Haar_DecompositionStep(vec, n, normal);

		for (int i = 0; i < n / 2; i++)
			_return[0][i + n / 2 - 1] = vec[i];

		for (int i = n / 2; i < n; i++)
			_return[1][i - 1] = vec[i];

		n /= 2;
	}

	return _return;
}

// Decomposição completa de tipo Haar de um vetor.
void Haar_Decomposition(double *vec, int n, bool normal)
{
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] / sqrt(float(n));
	
	while (n > 1)
	{
		Haar_DecompositionStep(vec, n, normal);
		n /= 2;
	}
}

void INT_Haar_Decomposition(interval *vec, int n, bool normal)
{
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] / sqrt(interval(n));
	
	while (n > 1)
	{
		INT_Haar_DecompositionStep(vec, n, normal);
		n /= 2;
	}
}

// Transformação Haar padrão de uma matriz.
void Haar_StandardDecomposition(double **matrix, int rows, int cols, bool normal)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		temp_row[j] = matrix[i][j];

		Haar_Decomposition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
		matrix[i][j] = temp_row[j];
	}
	
	for (int i = 0;i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		temp_col[j] = matrix[j][i];
			
		Haar_Decomposition(temp_col, rows, normal);
		
		for (int j = 0; j < rows; j++)
		matrix[j][i] = temp_col[j];
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

void INT_Haar_StandardDecomposition(interval **matrix, int rows, int cols, bool normal)
{
	interval *temp_row = new interval[cols];
	interval *temp_col = new interval[rows];
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		temp_row[j] = matrix[i][j];

		INT_Haar_Decomposition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
		matrix[i][j] = temp_row[j];
	}
	
	for (int i = 0;i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		temp_col[j] = matrix[j][i];
			
		INT_Haar_Decomposition(temp_col, rows, normal);
		
		for (int j = 0; j < rows; j++)
		matrix[j][i] = temp_col[j];
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

// Transformação Haar não-padrão de uma matriz.
void Haar_NonStandardDecomposition(double **matrix, int rows, int cols, bool normal)
{
	int h = rows, w = cols;
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];
	
	
	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j] / cols;
	}
	
	
	while (w > 1 || h > 1)
	{
		if (w > 1)
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				temp_row[j] = matrix[i][j];

			Haar_DecompositionStep(temp_row, w, normal);

			for (int j = 0; j < w; j++)
				matrix[i][j] = temp_row[j];
		}
		
		if (h > 1)
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
				temp_col[j] = matrix[j][i];
				
			Haar_DecompositionStep(temp_col, h, normal);
			
			for (int j = 0; j < h; j++)
				matrix[j][i] = temp_col[j];
		}
		
		if (w > 1) w /= 2;
		if (h > 1) h /= 2;
	}
}

void INT_Haar_NonStandardDecomposition(interval **matrix, int rows, int cols, bool normal)
{
	int h = rows, w = cols;
	interval *temp_row = new interval[cols];
	interval *temp_col = new interval[rows];
	
	
	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j] / interval(cols);
	}
	
	
	while (w > 1 || h > 1)
	{
		if (w > 1)
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				temp_row[j] = matrix[i][j];

			INT_Haar_DecompositionStep(temp_row, w, normal);

			for (int j = 0; j < w; j++)
				matrix[i][j] = temp_row[j];
		}
		
		if (h > 1)
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
				temp_col[j] = matrix[j][i];
				
			INT_Haar_DecompositionStep(temp_col, h, normal);
			
			for (int j = 0; j < h; j++)
				matrix[j][i] = temp_col[j];
		}
		
		if (w > 1) w /= 2;
		if (h > 1) h /= 2;
	}
}

// Executa o tratamento adequado para a matriz, baseado na HWT.
void Haar_MatrixDecomposition(double **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
		cout << "Não é uma matriz quadrada!" << endl;
		return;
	}
	*/
	
	if (standard) Haar_StandardDecomposition(matrix, rows, cols, normal);
	else          Haar_NonStandardDecomposition(matrix, rows, cols, normal);
}

void INT_Haar_MatrixDecomposition(interval **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
		cout << "Não é uma matriz quadrada!" << endl;
		return;
	}
	*/
	
	if (standard) INT_Haar_StandardDecomposition(matrix, rows, cols, normal);
	else          INT_Haar_NonStandardDecomposition(matrix, rows, cols, normal);
}

// Passo de composição tipo Haar de um vetor.
void Haar_CompositionStep(double *vec, int n, bool normal)
{
	double *vecp = new double[2*n];
	
	for (int i = 0; i < n; i++)
		vecp[i] = 0;
	
	for (int i = 0; i < n; i++)
	{
		vecp[2*i]     = vec[i] + vec[n + i];
		vecp[2*i + 1] = vec[i] - vec[n + i];
		
		if (normal)
		{
			vecp[2*i]     /= sqrt(2.0);
			vecp[2*i + 1] /= sqrt(2.0);
		}
	}
	
	for (int i = 0; i < 2*n; i++)
		vec[i] = vecp[i];
	
	delete [] vecp;
}

void INT_Haar_CompositionStep(interval *vec, int n, bool normal)
{
	interval *vecp = new interval[2*n];
	
	for (int i = 0; i < n; i++)
		vecp[i] = interval(0);
	
	for (int i = 0; i < n; i++)
	{
		vecp[2*i]     = vec[i] + vec[n + i];
		vecp[2*i + 1] = vec[i] - vec[n + i];
		
		if (normal)
		{
			vecp[2*i]     /= sqrt(interval(2));
			vecp[2*i + 1] /= sqrt(interval(2));
		}
	}
	
	for (int i = 0; i < 2*n; i++)
		vec[i] = vecp[i];
	
	delete [] vecp;
}

// Composição completa tipo Haar de um vetor.
void Haar_Composition(double *vec, int n, bool normal)
{
	for (int i = 1; i < n; i = i*2)
		Haar_CompositionStep(vec, i, normal);
	
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] * sqrt(float(n));
}

void INT_Haar_Composition(interval *vec, int n, bool normal)
{
	for (int i = 1; i < n; i = i*2)
		INT_Haar_CompositionStep(vec, i, normal);
	
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] * sqrt(interval(n));
}

// Transformação Inversa de Haar padrão de uma matriz.
void Haar_StandardComposition(double **matrix, int rows, int cols, bool normal)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];
	
	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		temp_col[j] = matrix[j][i];
			
		Haar_Composition(temp_col, rows, normal);
		
		for (int j = 0; j < rows; j++)
		matrix[j][i] = temp_col[j];
	}
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		temp_row[j] = matrix[i][j];

		Haar_Composition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
		matrix[i][j] = temp_row[j];
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

void INT_Haar_StandardComposition(interval **matrix, int rows, int cols, bool normal)
{
	interval *temp_row = new interval[cols];
	interval *temp_col = new interval[rows];
	
	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
		temp_col[j] = matrix[j][i];
			
		INT_Haar_Composition(temp_col, rows, normal);
		
		for (int j = 0; j < rows; j++)
		matrix[j][i] = temp_col[j];
	}
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		temp_row[j] = matrix[i][j];

		INT_Haar_Composition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
		matrix[i][j] = temp_row[j];
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

// Transformação Inversa de Haar não-padrão de uma matriz.
void Haar_NonStandardComposition(double **matrix, int rows, int cols, bool normal)
{
	int h = 1, w = 1;
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];
	
	while (w < cols || h < rows)
	{
		if (h < rows)
		for (int i = 0; i < 2*w; i++)
		{
			for (int j = 0; j < rows; j++)
				temp_col[j] = matrix[j][i];
			
			Haar_CompositionStep(temp_col, h, normal);
			
			for (int j = 0; j < rows; j++)
				matrix[j][i] = temp_col[j];
		}
		
		if (w < cols)
		for (int i = 0; i < 2*h; i++)
		{
			for (int j = 0; j < cols; j++)
				temp_row[j] = matrix[i][j];

			Haar_CompositionStep(temp_row, w, normal);

			for (int j = 0; j < cols; j++)
				matrix[i][j] = temp_row[j];
		}
		
		if (w < cols) w = w*2;
		if (h < rows) h = h*2;
	}
	
	
	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j]*cols;
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

void INT_Haar_NonStandardComposition(interval **matrix, int rows, int cols, bool normal)
{
	int h = 1, w = 1;
	interval *temp_row = new interval[cols];
	interval *temp_col = new interval[rows];
	
	while (w < cols || h < rows)
	{
		if (h < rows)
		for (int i = 0; i < 2*w; i++)
		{
			for (int j = 0; j < rows; j++)
				temp_col[j] = matrix[j][i];
			
			INT_Haar_CompositionStep(temp_col, h, normal);
			
			for (int j = 0; j < rows; j++)
				matrix[j][i] = temp_col[j];
		}
		
		if (w < cols)
		for (int i = 0; i < 2*h; i++)
		{
			for (int j = 0; j < cols; j++)
				temp_row[j] = matrix[i][j];

			INT_Haar_CompositionStep(temp_row, w, normal);

			for (int j = 0; j < cols; j++)
				matrix[i][j] = temp_row[j];
		}
		
		if (w < cols) w = w*2;
		if (h < rows) h = h*2;
	}
	
	
	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j]*cols;
	}
	
	delete [] temp_row;
	delete [] temp_col;
}

// Executa o tratamento adequado para a matriz, baseado na Inverse HWT.
void Haar_MatrixComposition(double **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
		cout << "Não é uma matriz quadrada!" << endl;
		return;
	}
	*/
	
	if (standard) Haar_StandardComposition(matrix, rows, cols, normal);
	else          Haar_NonStandardComposition(matrix, rows, cols, normal);
}

void INT_Haar_MatrixComposition(interval **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
		cout << "Não é uma matriz quadrada!" << endl;
		return;
	}
	*/
	
	if (standard) INT_Haar_StandardComposition(matrix, rows, cols, normal);
	else          INT_Haar_NonStandardComposition(matrix, rows, cols, normal);
}

real INT_diameter(interval x)
{
	real a1, a2;
	
	a1 = Inf(x);
	a2 = Sup(x);
	
	return a2 - a1;
}

real INT_error(interval *x, int n)
{
	real r = 0.0;
	
	for (int i = 0; i < n; i++)
	if (INT_diameter(x[i]) > r) r = INT_diameter(x[i]);
	
	return r;
}

real INT_error(interval **x, int linhas, int colunas)
{
	real r = 0.0, aux = 0.0;
	
	for (int i = 0; i < linhas; i++)
	{
		aux = INT_error(x[i], colunas);
		
		if (aux > r) r = aux;
	}
	
	return r;
}

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

    fscanf (imagem,"%s", auxc);
    if (auxc[0]=='#')
    {
        while(fgetc(imagem) != '\n');
        fscanf (imagem,"%d", &imgInfo.x);
    }
    else imgInfo.x = atoi(auxc);

	fscanf (imagem,"%d", &imgInfo.y);
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

void Haar_Compress(double **matrix, int rows, int cols, float percentege)
{
	double low, high, dif, cut, count;

	low = high = matrix[0][0];

	for (int i = 0; i < rows; i++)
	for (int j = 0; j < cols; j++)
	{
		if (matrix[i][j] < low) low = matrix[i][j];

		if (matrix[i][j] > high) high = matrix[i][j];
	}

	dif = high - low;

	if (dif == 0.0) return;

	cut = percentege * (dif / 100.0);

	cout << "high:\t" << high << endl;
	cout << "low:\t" << low << endl;
	cout << "dif:\t" << dif << endl;
	cout << "cut:\t" << cut << endl;

	count = 0.0;

	for (int i = 0; i < rows; i++)
	for (int j = 0; j < cols; j++)
	{
		if (abs(matrix[i][j]) <= cut)
		{
			count++;
			matrix[i][j] = 0.0;
		}
	}

	cout << "count:\t" << count << endl;
}