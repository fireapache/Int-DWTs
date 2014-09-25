#include "tests.h"

void matrixCopy(double **m1, double **m2, int x, int y)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			m2[i][j] = m1[i][j];
		}
	}
}

int test5(float percentage)
{
	int n = 8;

	double *dMat = new double[n];
	interval *iMat = new interval[n];

	for (int i = 0; i < n; ++i)
	{
		dMat[i] = i * i;
		iMat[i] = interval(i * i);
	}

	Haar_Decomposition(dMat, n, false);
	INT_Haar_Decomposition(iMat, n, false);

	printVector<double>(dMat, n);
	printVector<interval>(iMat, n);

	Haar_Compression(dMat, n, percentage);
	INT_Haar_Compression(iMat, n, percentage);

	printVector<double>(dMat, n);
	printVector<interval>(iMat, n);

	delete[] dMat;
	delete[] iMat;

	return 1;
}

int test6(float percentage)
{
	int n = 8;

	double **dMat = new double*[n];
	interval **iMat = new interval*[n];

	for (int i = 0; i < n; ++i)
	{
		dMat[i] = new double[n];
		iMat[i] = new interval[n];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			dMat[i][j] = double(i * j);
			iMat[i][j] = interval(i * j);
		}
	}

	printMatrix<double>(dMat, n);
	printMatrix<interval>(iMat, n);

	// Compression

	printMatrix<double>(dMat, n);
	printMatrix<interval>(iMat, n);

	for (int i = 0; i < n; ++i)
	{
		delete[] dMat[i];
		delete[] iMat[i];
	}

	delete[] dMat;
	delete[] iMat;

	return 1;
}

int test4(const char *ppmfilepath, double percentage)
{
	double **image = NULL;
	ImageInfo imageInfo;
	char *outppmpath;
	int charCount;

	charCount = 0;

	while (ppmfilepath[charCount] != '.')
		charCount++;

	if (charCount < 1) return 1;

	outppmpath = new char[charCount + 9];

	for (int i = 0; i < charCount ; ++i)
		outppmpath[i] = ppmfilepath[i];

	strcpy(&outppmpath[charCount], "_out.ppm");

	image = carregar_imagem((char*)(ppmfilepath), &imageInfo);

	Haar_MatrixDecomposition(image, imageInfo.x, imageInfo.y, true, false);

	Haar_PerLevel_Matrix_Compression(image, imageInfo.x, percentage);

	Haar_MatrixComposition(image, imageInfo.x, imageInfo.y, true, false);	

	escrever_imagem(outppmpath, image, imageInfo);

	for (int i = 0; i < imageInfo.x; i++)
	{
		delete[] image[i];
	}

	delete[] image;

	return 0;
}

int test3(const char *ppmfilepath, double percentage)
{
	double **image = NULL;
	ImageInfo imageInfo;
	char *outppmpath;
	int charCount;

	charCount = 0;

	while (ppmfilepath[charCount] != '.')
		charCount++;

	if (charCount < 1) return 1;

	outppmpath = new char[charCount + 9];

	for (int i = 0; i < charCount ; ++i)
		outppmpath[i] = ppmfilepath[i];

	strcpy(&outppmpath[charCount], "_out.ppm");

	image = carregar_imagem((char*)(ppmfilepath), &imageInfo);

	Haar_MatrixDecomposition(image, imageInfo.x, imageInfo.y, true, false);

	Haar_Matrix_Compression(image, imageInfo.x, percentage);

	Haar_MatrixComposition(image, imageInfo.x, imageInfo.y, true, false);	

	escrever_imagem(outppmpath, image, imageInfo);

	for (int i = 0; i < imageInfo.x; i++)
	{
		delete[] image[i];
	}

	delete[] image;

	return 0;
}

void test2_1()
{
	double **matrix = new double*[POINTS];
	double **test = new double*[POINTS];

	for (int i = 0; i < POINTS; i++)
	{
		matrix[i] = new double[POINTS];
		test[i] = new double[POINTS];
	}

	for (int i = 0; i < POINTS; i++)
	{
		for (int j = 0; j < POINTS; j++)
		{
			matrix[i][j] = (i + 1) + (j + 1);
			test[i][j] = (i + 1) + (j + 1);
		}
	}

	cout.setf(std::ios::fixed);
	cout << '\t';
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, false);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, false);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, false);

	cout << getTimeCounter() << '\t';

	cout << '\t';
	matrixCopy(test, matrix, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, true);

	cout << getTimeCounter() << '\t';
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, false);

	cout << getTimeCounter() << endl;
	cout.unsetf(std::ios::fixed);

	for (int i = 0; i < POINTS; i++)
	{
		delete[] matrix[i];
	}

	delete[] matrix;
}

void test2(unsigned int n)
{
	cout << endl << endl << "Decomposicao:";
	cout << '\t' << "MPN";
	cout << '\t' << "MPO";
	cout << '\t' << "MNPN";
	cout << '\t' << "MNPO";
	cout << '\t' << "Composicao";
	cout << '\t' << "MPN";
	cout << '\t' << "MPO";
	cout << '\t' << "MNPN";
	cout << '\t' << "MNPO";
	cout << endl << endl;

	for (unsigned int i = 0; i < n; i++)
	{
		test2_1();
	}
}

void test1()
{
	double **matrix = new double*[POINTS];
	double **test = new double*[POINTS];

	for (int i = 0; i < POINTS; i++)
	{
		matrix[i] = new double[POINTS];
		test[i] = new double[POINTS];
	}

	for (int i = 0; i < POINTS; i++)
	{
		for (int j = 0; j < POINTS; j++)
		{
			matrix[i][j] = (i + 1) + (j + 1);
			test[i][j] = (i + 1) + (j + 1);
		}
	}

	cout << endl << endl << "Decomposicao:" << endl << endl;
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, false);

	cout << "Metodo Padrao Novo: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, true);

	cout << "Metodo Padrao Original: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, false);

	cout << "Metodo Nao-Padrao Novo: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, false);

	cout << "Metodo Nao-Padrao Original: " << getTimeCounter() << endl;

	cout << endl << endl << "Composicao:" << endl << endl;
	matrixCopy(test, matrix, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, true);

	cout << "Metodo Padrao Novo: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, true);

	cout << "Metodo Padrao Original: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, true);

	cout << "Metodo Nao-Padrao Novo: " << getTimeCounter() << endl;
	matrixCopy(matrix, test, POINTS, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, false);

	cout << "Metodo Nao-Padrao Original: " << getTimeCounter() << endl;


	for (int i = 0; i < POINTS; i++)
	{
		delete[] matrix[i];
	}

	delete[] matrix;
}

void test0()
{
	double image[POINTS];
	double t1[POINTS], t2[POINTS], t3[POINTS];
	double x1 = 0.0, x2 = 1.0;
	double x, alpha;

	for (int i = 0; i < POINTS; i++)
	{
		alpha = (double)i / (double)(POINTS);
		x = x1 * (1.0 - alpha) + x2 * alpha;
		image[i] = sin(2 * PI * x) + 1 - sqrt(abs(8 * PI * (x - 0.5)));
		t1[i] = t2[i] = t3[i] = image[i];
	}

	gnuplot_dat_Vdecomposition("1_V.dat", x1, x2, image, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("1_1W.dat", x1, x2, image, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("1_2W.dat", x1, x2, image, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("1_3W.dat", x1, x2, image, POINTS, 3, true);

	Haar_Decomposition(t1, POINTS, true);
	//Haar_Compression(t1, POINTS, 0.01);
	Haar_PerLevel_Compression(t1, POINTS, 0.01);
	Haar_Composition(t1, POINTS, true);

	gnuplot_dat_Vdecomposition("2_V.dat", x1, x2, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_1W.dat", x1, x2, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_2W.dat", x1, x2, t1, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("2_3W.dat", x1, x2, t1, POINTS, 3, true);

	Haar_Decomposition(t2, POINTS, true);
	//Haar_Compression(t2, POINTS, 0.02);
	Haar_PerLevel_Compression(t2, POINTS, 0.02);
	Haar_Composition(t2, POINTS, true);

	gnuplot_dat_Vdecomposition("3_V.dat", x1, x2, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_1W.dat", x1, x2, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_2W.dat", x1, x2, t2, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("3_3W.dat", x1, x2, t2, POINTS, 3, true);

	Haar_Decomposition(t3, POINTS, true);
	//Haar_Compression(t3, POINTS, 0.035);
	Haar_PerLevel_Compression(t3, POINTS, 0.035);
	Haar_Composition(t3, POINTS, true);

	gnuplot_dat_Vdecomposition("4_V.dat", x1, x2, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_1W.dat", x1, x2, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_2W.dat", x1, x2, t3, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("4_3W.dat", x1, x2, t3, POINTS, 3, true);
}
