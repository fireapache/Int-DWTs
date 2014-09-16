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
	//Haar_Compress(t1, POINTS, 0.01);
	Haar_Levels_Compress(t1, POINTS, 0.01);
	Haar_Composition(t1, POINTS, true);

	gnuplot_dat_Vdecomposition("2_V.dat", x1, x2, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_1W.dat", x1, x2, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_2W.dat", x1, x2, t1, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("2_3W.dat", x1, x2, t1, POINTS, 3, true);

	Haar_Decomposition(t2, POINTS, true);
	//Haar_Compress(t2, POINTS, 0.02);
	Haar_Levels_Compress(t2, POINTS, 0.02);
	Haar_Composition(t2, POINTS, true);

	gnuplot_dat_Vdecomposition("3_V.dat", x1, x2, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_1W.dat", x1, x2, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_2W.dat", x1, x2, t2, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("3_3W.dat", x1, x2, t2, POINTS, 3, true);

	Haar_Decomposition(t3, POINTS, true);
	//Haar_Compress(t3, POINTS, 0.035);
	Haar_Levels_Compress(t3, POINTS, 0.035);
	Haar_Composition(t3, POINTS, true);

	gnuplot_dat_Vdecomposition("4_V.dat", x1, x2, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_1W.dat", x1, x2, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_2W.dat", x1, x2, t3, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("4_3W.dat", x1, x2, t3, POINTS, 3, true);
}
