#include "tests.h"

int test7(const char *filepath)
{
	double **oImg, **fImg = NULL;
	ImageInfo imageInfo;
	ImageQuality<double> imgQuality;

	oImg = carregar_imagem((char*)(filepath), &imageInfo);
	fImg = carregar_imagem((char*)(filepath), &imageInfo);

	Haar_MatrixDecomposition(fImg, imageInfo.x, imageInfo.y, true, false);
	Haar_MatrixComposition(fImg, imageInfo.x, imageInfo.y, true, false);

	imgQuality = imageQuality<double>(oImg, fImg, 255.0, imageInfo.x);

	cout << endl << endl;
	cout << "Algoritmos originais nao padrao: MSE, EUC, PSNR" << endl << endl;
	cout << imgQuality.mse << endl << imgQuality.euc << endl << imgQuality.psnr << endl << endl;

	deleteMatrix<double>(fImg, imageInfo.x);
	fImg = carregar_imagem((char*)(filepath), &imageInfo);

	Haar_MatrixDecomposition(fImg, imageInfo.x, imageInfo.y, false, false);
	VinisMatrixNormalization(fImg, imageInfo.x, false);
	VinisMatrixNormalization(fImg, imageInfo.x, false, true);
	Haar_MatrixComposition(fImg, imageInfo.x, imageInfo.y, false, false);

	imgQuality = imageQuality<double>(oImg, fImg, 255.0, imageInfo.x);

	cout << "Algoritmos novos nao padrao: MSE, EUC, PSNR" << endl << endl;
	cout << imgQuality.mse << endl << imgQuality.euc << endl << imgQuality.psnr << endl << endl;

	deleteMatrix<double>(fImg, imageInfo.x);
	fImg = carregar_imagem((char*)(filepath), &imageInfo);

	Haar_MatrixDecomposition(fImg, imageInfo.x, imageInfo.y, true, true);
	Haar_MatrixComposition(fImg, imageInfo.x, imageInfo.y, true, true);

	imgQuality = imageQuality<double>(oImg, fImg, 255.0, imageInfo.x);

	cout << "Algoritmos originais padrao: MSE, EUC, PSNR" << endl << endl;
	cout << imgQuality.mse << endl << imgQuality.euc << endl << imgQuality.psnr << endl << endl;

	deleteMatrix<double>(fImg, imageInfo.x);
	fImg = carregar_imagem((char*)(filepath), &imageInfo);

	Haar_MatrixDecomposition(fImg, imageInfo.x, imageInfo.y, false, true);
	VinisMatrixNormalization(fImg, imageInfo.x, true);
	VinisMatrixNormalization(fImg, imageInfo.x, true, true);
	Haar_MatrixComposition(fImg, imageInfo.x, imageInfo.y, false, true);

	imgQuality = imageQuality<double>(oImg, fImg, 255.0, imageInfo.x);

	cout << "Algoritmos novos padrao: MSE, EUC, PSNR" << endl << endl;
	cout << imgQuality.mse << endl << imgQuality.euc << endl << imgQuality.psnr << endl << endl;

	deleteMatrix<double>(fImg, imageInfo.x);
	deleteMatrix<double>(oImg, imageInfo.x);

	return 0;
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

	Haar_MatrixDecomposition(dMat, n, n, false, false);
	INT_Haar_MatrixDecomposition(iMat, n, n, false, false);

	printMatrix<double>(dMat, n);
	printMatrix<interval>(iMat, n);

	Haar_Matrix_Compression(dMat, n, percentage);
	INT_Haar_Matrix_Compression(iMat, n, percentage);

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
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, false);

	cout << getTimeCounter() << '\t';
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, false);

	cout << getTimeCounter() << '\t';

	cout << '\t';
	copyMatrix(test, matrix, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, true);

	cout << getTimeCounter() << '\t';
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, true);

	cout << getTimeCounter() << '\t';
	copyMatrix(matrix, test, POINTS);
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
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, true);

	cout << "Metodo Padrao Original: " << getTimeCounter() << endl;
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, false);

	cout << "Metodo Nao-Padrao Novo: " << getTimeCounter() << endl;
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixDecomposition(test, POINTS, POINTS, true, false);

	cout << "Metodo Nao-Padrao Original: " << getTimeCounter() << endl;

	cout << endl << endl << "Composicao:" << endl << endl;
	copyMatrix(test, matrix, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, true);
	VinisMatrixNormalization(test, POINTS, true, true);

	cout << "Metodo Padrao Novo: " << getTimeCounter() << endl;
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, true, true);

	cout << "Metodo Padrao Original: " << getTimeCounter() << endl;
	copyMatrix(matrix, test, POINTS);
	startTimeCounter();

	Haar_MatrixComposition(test, POINTS, POINTS, false, false);
	VinisMatrixNormalization(test, POINTS, false, true);

	cout << "Metodo Nao-Padrao Novo: " << getTimeCounter() << endl;
	copyMatrix(matrix, test, POINTS);
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

void fundamentalTest0()
{
	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Entrada: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;
	
	double **mat = new double*[4];
	
	for (int m = 0; m < 4; m++)
		mat[m] = new double[4];
	
	mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2;
	mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5;
	mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1;
	mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1;

	interval **intmat = new interval*[4];
	
	for (int m = 0; m < 4; m++)
		intmat[m] = new interval[4];
	
	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	{
		intmat[i][j] = interval(mat[i][j]);
	}
	
	double **matriz = new double*[4];
	
	for (int m = 0; m < 4; m++)
		matriz[m] = new double[4];
	
	copyMatrix<double>(mat, matriz, 4);

	interval **intmatriz = new interval*[4];
	
	for (int m = 0; m < 4; m++)
		intmatriz[m] = new interval[4];
	
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	cout << "==========================================" << endl;
	cout << "(Algoritmos originais) Processo não normalizado: \n" << endl;
	
	cout << "Pontual padrão: " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, false, true);
	printMatrix<double>(matriz, 4);
	Haar_MatrixComposition(matriz, 4, 4, false, true);
	printMatrix<double>(matriz, 4);
	copyMatrix<double>(mat, matriz, 4);
	
	cout << "Pontual não-padrão: " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, false, false);
	printMatrix<double>(matriz, 4);
	Haar_MatrixComposition(matriz, 4, 4, false, false);
	printMatrix<double>(matriz, 4);
	cout << endl;
	
	cout << "Intervalar padrão: " << endl << endl;
	
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, false, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_Haar_MatrixComposition(intmatriz, 4, 4, false, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	
	cout << endl;
	
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	cout << "Intervalar não-padrão: " << endl << endl;
	
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, false, false);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_Haar_MatrixComposition(intmatriz, 4, 4, false, false);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "(Algoritmos originais) Processo normalizado: \n" << endl;
	
	copyMatrix<double>(mat, matriz, 4);
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	cout << "Pontual padrão (normalizado): " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, true, true);
	printMatrix<double>(matriz, 4);
	Haar_MatrixComposition(matriz, 4, 4, true, true);
	printMatrix<double>(matriz, 4);
	
	cout << endl;
	
	copyMatrix<double>(mat, matriz, 4);
	
	cout << "Pontual não-padrão (normalizado): " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, true, false);
	printMatrix<double>(matriz, 4);
	Haar_MatrixComposition(matriz, 4, 4, true, false);
	printMatrix<double>(matriz, 4);
	
	cout << endl;
	
	cout << "Intervalar padrão (normalizado): " << endl << endl;
	
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, true, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_Haar_MatrixComposition(intmatriz, 4, 4, true, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	
	cout << endl;
	
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	cout << "Intervalar não-padrão (normalizado): " << endl << endl;
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, true, false);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_Haar_MatrixComposition(intmatriz, 4, 4, true, false);
	printMatrix<interval>(intmatriz, 4);
	
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;

	cout << endl;

	cout << "==========================================" << endl;
	cout << "(Algoritmos novos) Processo normalizado: \n" << endl;

	copyMatrix<double>(mat, matriz, 4);
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	cout << "Pontual padrão (normalizado): " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, false, true);
	VinisMatrixNormalization(matriz, 4, true);
	printMatrix<double>(matriz, 4);
	VinisMatrixNormalization(matriz, 4, true, true);
	Haar_MatrixComposition(matriz, 4, 4, false, true);
	printMatrix<double>(matriz, 4);
	
	cout << endl;
	
	copyMatrix<double>(mat, matriz, 4);
	
	cout << "Pontual não-padrão (normalizado): " << endl << endl;

	Haar_MatrixDecomposition(matriz, 4, 4, false, false);
	VinisMatrixNormalization(matriz, 4, false);
	printMatrix<double>(matriz, 4);
	VinisMatrixNormalization(matriz, 4, false, true);
	Haar_MatrixComposition(matriz, 4, 4, false, false);
	printMatrix<double>(matriz, 4);
	
	cout << endl;
	
	cout << "Intervalar padrão (normalizado): " << endl << endl;
	
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, false, true);
	INT_VinisMatrixNormalization(intmatriz, 4, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_VinisMatrixNormalization(intmatriz, 4, true, true);
	INT_Haar_MatrixComposition(intmatriz, 4, 4, false, true);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	
	cout << endl;
	
	copyMatrix<interval>(intmat, intmatriz, 4);
	
	INT_Haar_MatrixDecomposition(intmatriz, 4, 4, false, false);
	INT_VinisMatrixNormalization(intmatriz, 4, false);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	INT_VinisMatrixNormalization(intmatriz, 4, false, true);
	INT_Haar_MatrixComposition(intmatriz, 4, 4, false, false);
	printMatrix<interval>(intmatriz, 4);
	cout << '\t' << "Intervalo de erro:  " << INT_error(intmatriz, 4, 4) << endl << endl;
	
	cout << endl;
}

void fundamentalTest1()
{
	double *vec = new double[8];
	double *vector = new double[8];
	double **resultBuffer = new double*[8];
	interval *intVec = new interval[8];
	interval *intVector = new interval[8];
	interval **intResultBuffer = new interval*[4];

	vector[0] = 1.0;  vector[1] = 4.0;  vector[2] = 5.0; vector[3] = 2.0;
	vector[4] = 10.0; vector[5] = -3.0; vector[6] = 7.0; vector[7] = 6.0;

	for (int i = 0; i < 8; i++) intVector[i] = interval(vector[i]);

	for (int i = 0; i < 4; ++i)
	{
		resultBuffer[i] = new double[8];
		intResultBuffer[i] = new interval[8];
	}

	copyVector<double>(vector, vec, 8);
	copyVector<interval>(intVector, intVec, 8);

	copyVector<double>(vec, resultBuffer[0], 8);
	copyVector<interval>(intVec, intResultBuffer[0], 8);

	Haar_Decomposition(vec, 8, true);
	INT_Haar_Decomposition(intVec, 8, true);

	copyVector<double>(vec, resultBuffer[1], 8);
	copyVector<interval>(intVec, intResultBuffer[1], 8);

	Haar_Compression(vec, 8, 0.01);
	INT_Haar_Compression(intVec, 8, 0.01);

	copyVector<double>(vec, resultBuffer[2], 8);
	copyVector<interval>(intVec, intResultBuffer[2], 8);

	Haar_Composition(vec, 8, true);
	INT_Haar_Composition(intVec, 8, true);

	copyVector<double>(vec, resultBuffer[3], 8);
	copyVector<interval>(intVec, intResultBuffer[3], 8);

	cout << endl;
	cout << "========= Fundamental Test 1 =========" << endl;
	cout << endl;
	cout << "========= Testing compress methods on original and new algorithms ==========" << endl;
	cout << endl;
	cout << "========= Using original algorithms:";
	cout << endl;
	cout << endl;	

	printVectors<double>(resultBuffer, 8, 4, 10);
	cout << endl;
	printVectors<interval>(intResultBuffer, 8, 4, 1);

	copyVector<double>(vector, vec, 8);
	copyVector<interval>(intVector, intVec, 8);

	copyVector<double>(vec, resultBuffer[0], 8);
	copyVector<interval>(intVec, intResultBuffer[0], 8);

	Haar_Decomposition(vec, 8, false);
	VinisNormalization(vec, 8);
	INT_Haar_Decomposition(intVec, 8, false);
	INT_VinisNormalization(intVec, 8);

	copyVector<double>(vec, resultBuffer[1], 8);
	copyVector<interval>(intVec, intResultBuffer[1], 8);

	Haar_Compression(vec, 8, 0.01);
	INT_Haar_Compression(intVec, 8, 0.01);

	copyVector<double>(vec, resultBuffer[2], 8);
	copyVector<interval>(intVec, intResultBuffer[2], 8);

	VinisNormalization(vec, 8, true);
	Haar_Composition(vec, 8, false);
	INT_VinisNormalization(intVec, 8, true);
	INT_Haar_Composition(intVec, 8, false);

	copyVector<double>(vec, resultBuffer[3], 8);
	copyVector<interval>(intVec, intResultBuffer[3], 8);

	cout << endl;
	cout << "========= Using new algorithms:";
	cout << endl;
	cout << endl;

	printVectors<double>(resultBuffer, 8, 4, 10);
	cout << endl;
	printVectors<interval>(intResultBuffer, 8, 4, 1);
	
}

void fundamentalTest(unsigned int n)
{
	switch(n)
	{
		case 0:
			fundamentalTest0();
			break;
		case 1:
			fundamentalTest1();
			break;
		default:
			cout << "\n\tFundamental test " << n << " not found!\n\n";
	}
}
