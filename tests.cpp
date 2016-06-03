#include "tests.h"

void test15Desc()
{
	cout << endl;
	cout << "	==================== Test 4 ====================" << endl;
	cout << endl;
	cout << "	* This test is about mesurement of time and calculus " << endl;
	cout << "	exactitude for the two-dimensional À-Trous Haar Wavelet " << endl;
	cout << "	Transform." << endl;
}

void test15Param()
{
	cout << endl;
	cout << "	Parameters: <order> <levels>" << endl;
	cout << "	" << endl;
	cout << "	<order> must be an unsigned integer greater " << endl;
	cout << "	than 0 and also be power of two." << endl;
	cout << "	" << endl;
	cout << "	<levels> must be an unsigned integet and " << endl;
	cout << "	be also greater than 0." << endl;
	cout << endl;
}

real test15_ProcessError(interval **imat, int n, uint levels, uint opt, bool standard, bool old)
{
	interval ***data = NULL;
	interval **comp = NULL;
	real error = real(0);
#ifndef WIN32
	if (opt == 0)
	{
		data = Haar_atrous_MatrixDecomposition(imat, n, n, levels, old, standard);
		if (!old) Haar_atrous_MatrixNormalization(imat, data, n, n, levels);
	}
	else if (opt == 1)
	{
		data = genRandomMatrices<interval>(levels * 4, n, n, n);
		comp = Haar_atrous_MatrixComposition(data, n, n, levels, standard);
	}
	else
	{
		data = Haar_atrous_MatrixDecomposition(imat, n, n, levels, old, standard);
		if (!old) Haar_atrous_MatrixNormalization(imat, data, n, n, levels);
		comp = Haar_atrous_MatrixComposition(data, n, n, levels, standard);
	}

	if (opt == 0)
	{
		error = INT_error(data, levels * 4, n, n);
		deleteMatrices(data, levels * 4, n);
	}
	else
	{
		error = INT_error(comp, n, n);
		deleteMatrices(data, levels * 4, n);
		deleteMatrix(comp, n);
	}
#endif
	return error;
}

double test15_Process(double **mat, int n, uint levels, uint opt, bool standard, bool old, double **&outMat)
{
	double ***data = NULL;
	double **comp = NULL;
	double time;

	startTimeCounter();

	if (opt == 0)
	{
		data = Haar_atrous_MatrixDecomposition(mat, n, n, levels, old, standard);
		if (!old) Haar_atrous_MatrixNormalization(mat, data, n, n, levels);
	}
	else if (opt == 1)
	{
		data = genRandomMatrices<double>(levels * 4, n, n, n);
		comp = Haar_atrous_MatrixComposition(data, n, n, levels, standard);
	}
	else if (opt == 2)
	{
		data = Haar_atrous_MatrixDecomposition(mat, n, n, levels, old, standard);
		if (!old) Haar_atrous_MatrixNormalization(mat, data, n, n, levels);
		comp = Haar_atrous_MatrixComposition(data, n, n, levels, standard);
	}

	time = getTimeCounter();
	outMat = comp;

	if (data != NULL) deleteMatrices(data, levels * 4, n);
	if (opt != 2 && comp != NULL) deleteMatrix(comp, n);

	return time;
}

void test15_Script(double **inputMat, interval **inputiMat, uint n, uint levels, bool standard)
{
	double **comp;
	double *times = new double[30];			// Vector to store the time of each execution.
	real error, oldError;
	TimeMesurement timeMesure, oldMesure;

	if (standard) cout << "========== Standard Decomposition: " << endl;
	else cout << "========== Non Standard Decomposition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;

	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 0, standard, true, comp);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	oldError = test15_ProcessError(inputiMat, n, levels, 0, standard, true);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 0, standard, false, comp);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	error = test15_ProcessError(inputiMat, n, levels, 0, standard, false);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	cout << endl;
	if (standard) cout << "========== Standard Composition: " << endl;
	else cout << "========== Non Standard Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << endl;
	
	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 1, standard, true, comp);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	oldError = test15_ProcessError(inputiMat, n, levels, 1, standard, true);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 1, standard, false, comp);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	error = test15_ProcessError(inputiMat, n, levels, 1, standard, false);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	ImageQuality<double> imgQ, oldImgQ;

	cout << endl;
	if (standard) cout << "========== Standard Decomposition & Composition: " << endl;
	else cout << "========== Non Standard Decomposition & Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "EUC" << "\t\t";
	cout << "MSE" << "\t\t";
	cout << "PSNR" << endl;

	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 2, standard, true, comp);
		deleteMatrix(comp, n);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	oldError = test15_ProcessError(inputiMat, n, levels, 2, standard, true);

	cout << '\t' << oldError;
	
	test15_Process(inputMat, n, levels, 2, standard, true, comp);

	double oldMax = maxValue(comp, n, n);

	oldImgQ = imgQ = imageQuality(comp, inputMat, oldMax, n);

	delete [] comp;

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	for (int i = 0; i < 30; ++i)
	{
		times[i] = test15_Process(inputMat, n, levels, 2, standard, false, comp);
		deleteMatrix(comp, n);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	error = test15_ProcessError(inputiMat, n, levels, 2, standard, false);

	cout << '\t' << error;
	
	test15_Process(inputMat, n, levels, 2, standard, false, comp);

	double max = maxValue(comp, n, n);

	imgQ = imageQuality(comp, inputMat, max, n);
	
	delete [] comp;

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	cout << endl;
	cout << "Quality gain (%)" << endl;
	cout << endl;

	cout << "EUC" << "\t\t" << relativeGain(imgQ.euc, oldImgQ.euc) << endl;
	cout << "MSE" << "\t\t" << relativeGain(imgQ.mse, oldImgQ.mse) << endl;
	cout << "PSNR" << "\t\t" << -relativeGain(imgQ.psnr, oldImgQ.psnr) << endl;
	cout << "Error" << "\t\t" << relativeGain(error, oldError) << endl;

	cout << endl;

	delete [] times;
}

int test15(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << endl;
		cout << "	ERROR: lack of parameters." << endl;
		test15Desc();
		test15Param();
		return 1;
	}

	uint n = atoi(argv[0]);
	uint levels = atoi(argv[1]);

	if (n <= 0 || levels <= 0)
	{
		cout << endl;
		cout << "	ERROR: <order> and <levels> have to be greater than 0." << endl;
		test15Desc();
		test15Param();
		return 1;
	}
	else if (!isPowerOfTwo(n))
	{
		cout << endl;
		cout << "	ERROR: <order> has to be power of two." << endl;
		test15Desc();
		test15Param();
		return 1;
	}

	test15Desc();
	cout << endl;
	cout << "Computing " << levels << "levels of a matrix of order " << n  << "..." << endl;
	cout << n << " X " << n << " = " << n * n << "..." << endl;
	cout << endl;

	double **inputMat = new double*[n];		// Matrices used as input for all executions.
	interval **inputiMat = new interval*[n];

	double **auxMat = new double*[n];		// Auxiliary matrices used for each execution.
	interval **auxiMat = new interval*[n];

	// Allocate and fill input matrices with n*n random values.
	for (uint i = 0; i < n; ++i)
	{
		inputMat[i] = new double[n];
		inputiMat[i] = new interval[n];

		auxMat[i] = new double[n];
		auxiMat[i] = new interval[n];

		for (uint j = 0; j < n; ++j)
		{
			inputMat[i][j] = rand() % n;
			inputiMat[i][j] = interval(inputMat[i][j]);
		}
	}

	test15_Script(inputMat, inputiMat, n, levels, true);
	test15_Script(inputMat, inputiMat, n, levels, false);

	// Deallocating memory.

	for (uint i = 0; i < n; ++i)
	{
		delete [] inputMat[n];
		delete [] inputiMat[n];
		delete [] auxMat[n];
		delete [] auxiMat[n];
	}

	delete [] inputMat;
	delete [] inputiMat;
	delete [] auxMat;
	delete [] auxiMat;

	return 0;
}

void test14Desc()
{
	cout << endl;
	cout << "	==================== Test 3 ====================" << endl;
	cout << endl;
	cout << "	* This test is about mesurement of time and calculus " << endl;
	cout << "	exactitude for the one-dimensional À-Trous Haar Wavelet " << endl;
	cout << "	Transform." << endl;
}

void test14Param()
{
	cout << endl;
	cout << "	Parameters: <size> <levels>" << endl;
	cout << "	" << endl;
	cout << "	<size> must be an unsigned integer greater " << endl;
	cout << "	than 0 and also be power of two." << endl;
	cout << "	" << endl;
	cout << "	<levels> must be an unsigned integet and " << endl;
	cout << "	be also greater than 0." << endl;
	cout << endl;
}

real test14_OriginalDecompError(interval *ivec, int n, uint levels)
{
	interval **data;
	real error;

	data = Haar_atrous_Decomposition(ivec, n, levels, true);
	error = INT_error(data, levels, n);
	deleteMatrix(data, levels * 2);

	return error;
}

double test14_OriginalDecomp(double *vec, int n, uint levels)
{
	double timer = 0.0;
	double **data;

	// Get time of execution.
	startTimeCounter();
	data = Haar_atrous_Decomposition(vec, n, levels, true);
	timer = getTimeCounter();
	
	deleteMatrix(data, levels * 2);

	return timer;
}

real test14_DevelopedDecompError(interval *ivec, int n, uint levels)
{
	interval **data;
	real error;

	data = Haar_atrous_Decomposition(ivec, n, levels, false);
	Haar_atrous_Normalization(ivec, data, n, levels);
	error = INT_error(data, levels, n);
	deleteMatrix(data, levels * 2);

	return error;
}

double test14_DevelopedDecomp(double *vec, int n, uint levels)
{
	double timer = 0.0;
	double **data;

	// Get time of execution.
	startTimeCounter();
	data = Haar_atrous_Decomposition(vec, n, levels, false);
	Haar_atrous_Normalization(vec, data, n, levels);
	timer = getTimeCounter();
	
	deleteMatrix(data, levels * 2);

	return timer;
}

real test14_OriginalCompError(interval *ivec, uint n, uint levels)
{
	interval **data;
	interval *comp;
	real result;

	// Get time of execution.
	data = genRandomMatrix<interval>(levels * 2, n, n);
	comp = Haar_atrous_Composition(data, n, levels);
	result = INT_error(comp, n);

	deleteMatrix(data, levels * 2);
	delete [] comp;

	return result;
}

double test14_OriginalComp(double *vec, int n, uint levels)
{
	double timer = 0.0;
	double **data;
	double *result;

	// Get time of execution.
	data = genRandomMatrix<double>(levels * 2, n, n);
	startTimeCounter();
	result = Haar_atrous_Composition(data, n, levels);
	timer = getTimeCounter();

	deleteMatrix(data, levels * 2);
	delete [] result;

	return timer;
}

real test14_DevelopedCompError(interval *ivec, uint n, uint levels)
{
	interval **data;
	interval *comp;
	real result;

	// Get time of execution.
	data = genRandomMatrix<interval>(levels * 2, n, n);
	comp = Haar_atrous_Composition(data, n, levels);
	result = INT_error(comp, n);

	deleteMatrix(data, levels * 2);
	delete [] comp;

	return result;
}

double test14_DevelopedComp(double *vec, int n, uint levels)
{
	double timer = 0.0;
	double **data;
	double *result;

	// Get time of execution.
	data = genRandomMatrix<double>(levels * 2, n, n);
	startTimeCounter();
	result = Haar_atrous_Composition(data, n, levels);
	timer = getTimeCounter();

	deleteMatrix(data, levels * 2);
	delete [] result;

	return timer;
}

real test14_OriginalDecompCompError(interval *ivec, int n, uint levels)
{
	interval **data;
	interval *comp;
	real result;

	// Get time of execution.
	data = Haar_atrous_Decomposition(ivec, n, levels, true);
	comp = Haar_atrous_Composition(data, n, levels);
	result = INT_error(comp, n);

	deleteMatrix(data, levels * 2);
	delete [] comp;

	return result;
}

double test14_OriginalDecompComp(double *vec, double *&outVec, int n, uint levels, bool out = false)
{
	double timer = 0.0;
	double **data;
	double *result;

	// Get time of execution.
	startTimeCounter();
	data = Haar_atrous_Decomposition(vec, n, levels, true);
	result = Haar_atrous_Composition(data, n, levels);
	timer = getTimeCounter();

	deleteMatrix(data, levels * 2);

	if (out) outVec = result;
	else delete [] result;

	return timer;
}

real test14_DevelopedDecompCompError(interval *ivec, int n, uint levels)
{
	interval **data;
	interval *comp;
	real result;

	// Get time of execution.
	data = Haar_atrous_Decomposition(ivec, n, levels, false);
	Haar_atrous_Normalization(ivec, data, n, levels);
	comp = Haar_atrous_Composition(data, n, levels);

	result = INT_error(comp, n);

	deleteMatrix(data, levels * 2);
	delete [] comp;

	return result;
}

double test14_DevelopedDecompComp(double *vec, double *&outVec, int n, uint levels, bool out = false)
{
	double timer = 0.0;
	double **data;
	double *result;

	// Get time of execution.
	startTimeCounter();
	data = Haar_atrous_Decomposition(vec, n, levels, false);
	Haar_atrous_Normalization(vec, data, n, levels);
	result = Haar_atrous_Composition(data, n, levels);
	timer = getTimeCounter();

	deleteMatrix(data, levels * 2);

	if (out) outVec = result;
	else delete [] result;

	return timer;
}

int test14(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << endl;
		cout << "	ERROR: lack of parameters." << endl;
		test14Desc();
		test14Param();
		return 1;
	}

	uint n = atoi(argv[0]);
	uint levels = atoi(argv[1]);

	if (n <= 0 || levels <= 0)
	{
		cout << endl;
		cout << "	ERROR: <size> and <levels> have to be greater than 0." << endl;
		test14Desc();
		test14Param();
		return 1;
	}
	else if (!isPowerOfTwo(n))
	{
		cout << endl;
		cout << "	ERROR: <size> has to be power of two." << endl;
		test14Desc();
		test14Param();
		return 1;
	}

	test14Desc();
	cout << endl;
	cout << "Computing " << levels << " levels for a vector of size " << n  << "..." << endl;
	cout << endl;

	TimeMesurement timeMesure, oldMesure;

	double *inputVec = new double[n];		// Vectors used as input for all executions.
	interval *inputiVec = new interval[n];

	double *auxVec = new double[n];			// Auxiliary vectors used for each execution.
	interval *auxiVec = new interval[n];

	double *times = new double[30];			// Vector to store the time of each execution.

	real error, oldError;

	// Fill input vectors with n random values.
	for (uint i = 0; i < n; ++i)
	{
		inputVec[i] = rand() %  n;
		inputiVec[i] = interval(inputVec[i]);
	}

	cout << "========== Decomposition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_OriginalDecomp(auxVec, n, levels);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	oldError = test14_OriginalDecompError(auxiVec, n, levels);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_DevelopedDecomp(auxVec, n, levels);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	error = test14_DevelopedDecompError(auxiVec, n, levels);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	cout << endl;
	cout << "========== Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;
	
	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_OriginalComp(auxVec, n, levels);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	oldError = test14_OriginalCompError(auxiVec, n, levels);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_DevelopedComp(auxVec, n, levels);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	error = test14_DevelopedCompError(auxiVec, n, levels);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	ImageQuality<double> imgQ, oldImgQ;
	double *result;

	cout << endl;
	cout << "========== Decomposition & Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "EUC" << "\t\t";
	cout << "MSE" << "\t\t";
	cout << "PSNR" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_OriginalDecompComp(auxVec, result, n, levels);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	oldError = test14_OriginalDecompCompError(auxiVec, n, levels);

	cout << '\t' << oldError;
	
	test14_OriginalDecompComp(inputVec, result, n, levels, true);

	double oldMax = maxValue(result, n);

	oldImgQ = imgQ = imageQuality(result, inputVec, oldMax, n);

	delete [] result;

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test14_DevelopedDecompComp(auxVec, result, n, levels);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	error = test14_DevelopedDecompCompError(auxiVec, n, levels);

	cout << '\t' << error;
	
	test14_DevelopedDecompComp(inputVec, result, n, levels, true);

	double max = maxValue(result, n);

	imgQ = imageQuality(result, inputVec, max, n);

	delete [] result;

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	cout << endl;
	cout << "Quality gain (%)" << endl;
	cout << endl;

	cout << "EUC" << "\t\t" << relativeGain(imgQ.euc, oldImgQ.euc) << endl;
	cout << "MSE" << "\t\t" << relativeGain(imgQ.mse, oldImgQ.mse) << endl;
	cout << "PSNR" << "\t\t" << -relativeGain(imgQ.psnr, oldImgQ.psnr) << endl;
	cout << "Error" << "\t\t" << relativeGain(error, oldError) << endl;

	cout << endl;

	// Deallocating memory.
	delete [] inputVec;
	delete [] inputiVec;
	delete [] auxVec;
	delete [] auxiVec;
	delete [] times;

	return 0;
}

void test13Desc()
{
	cout << endl;
	cout << "	==================== Test 2 ====================" << endl;
	cout << endl;
	cout << "	* This test is about mesurement of time and calculus " << endl;
	cout << "	exactitude for the two-dimensional Haar Wavelet Transform." << endl;
}

void test13Param()
{
	cout << endl;
	cout << "	Parameters: <order>" << endl;
	cout << "	" << endl;
	cout << "	The <order> of matrix must be an unsigned integer " << endl;
	cout << "	greater than 0 and also be power of two." << endl;
	cout << endl;
}

double test13_OriginalDecomp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_MatrixDecomposition(mat, n, n, true, standard);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_MatrixDecomposition(imat, n, n, true, standard);
	}	

	return timer;
}

double test13_DevelopedDecomp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_MatrixDecomposition(mat, n, n, false, standard);
		if (standard) VinisStandardMatrixNormalization(mat, n);
		else VinisNonStandardMatrixNormalization(mat, n);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_MatrixDecomposition(imat, n, n, false, standard);
		if (standard) INT_VinisStandardMatrixNormalization(imat, n);
		else INT_VinisNonStandardMatrixNormalization(imat, n);
	}	

	return timer;
}

double test13_OriginalComp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_MatrixComposition(mat, n, n, true, standard);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_MatrixComposition(imat, n, n, true, standard);
	}	

	return timer;
}

double test13_DevelopedComp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		if (standard) VinisStandardMatrixNormalization(mat, n, true);
		else VinisNonStandardMatrixNormalization(mat, n, true);
		Haar_MatrixComposition(mat, n, n, false, standard);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		if (standard) INT_VinisStandardMatrixNormalization(imat, n, true);
		else INT_VinisNonStandardMatrixNormalization(imat, n, true);
		INT_Haar_MatrixComposition(imat, n, n, false, standard);
	}	

	return timer;
}

double test13_OriginalDecompComp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_MatrixDecomposition(mat, n, n, true, standard);
		Haar_MatrixComposition(mat, n, n, true, standard);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_MatrixDecomposition(imat, n, n, true, standard);
		INT_Haar_MatrixComposition(imat, n, n, true, standard);
	}	

	return timer;
}

double test13_DevelopedDecompComp(double **mat, interval **imat, int n, bool standard, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_MatrixDecomposition(mat, n, n, false, standard);

		if (standard)
		{
			VinisStandardMatrixNormalization(mat, n);
			VinisStandardMatrixNormalization(mat, n, true);
		}
		else
		{
			VinisNonStandardMatrixNormalization(mat, n);
			VinisNonStandardMatrixNormalization(mat, n, true);
		}

		Haar_MatrixComposition(mat, n, n, false, standard);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_MatrixDecomposition(imat, n, n, false, standard);
		
		if (standard)
		{
			INT_VinisStandardMatrixNormalization(imat, n);
			INT_VinisStandardMatrixNormalization(imat, n, true);
		}
		else
		{
			INT_VinisNonStandardMatrixNormalization(imat, n);
			INT_VinisNonStandardMatrixNormalization(imat, n, true);
		}

		INT_Haar_MatrixComposition(imat, n, n, false, standard);
	}	

	return timer;
}

void test13_Script(double **inputMat, double **auxMat, interval **inputiMat, interval **auxiMat, uint n, bool standard)
{
	double *times = new double[30];			// Vector to store the time of each execution.
	TimeMesurement timeMesure, oldMesure;
	real error, oldError;

	if (standard) cout << "========== Standard Decomposition: " << endl;
	else cout << "========== Non Standard Decomposition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_OriginalDecomp(auxMat, auxiMat, n, standard, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_OriginalDecomp(auxMat, auxiMat, n, standard, false);

	oldError = INT_error(auxiMat, n, n);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_DevelopedDecomp(auxMat, auxiMat, n, standard, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_DevelopedDecomp(auxMat, auxiMat, n, standard, false);

	error = INT_error(auxiMat, n, n);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	cout << endl;
	if (standard) cout << "========== Standard Composition: " << endl;
	else cout << "========== Non Standard Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;
	
	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_OriginalComp(auxMat, auxiMat, n, standard, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_OriginalComp(auxMat, auxiMat, n, standard, false);

	oldError = INT_error(auxiMat, n, n);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_DevelopedComp(auxMat, auxiMat, n, standard, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_DevelopedComp(auxMat, auxiMat, n, standard, false);

	error = INT_error(auxiMat, n, n);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	ImageQuality<double> imgQ, oldImgQ;

	cout << endl;
	if (standard) cout << "========== Standard Decomposition & Composition: " << endl;
	else cout << "========== Non Standard Decomposition & Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "EUC" << "\t\t";
	cout << "MSE" << "\t\t";
	cout << "PSNR" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_OriginalDecompComp(auxMat, auxiMat, n, standard, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_OriginalDecompComp(auxMat, auxiMat, n, standard, false);

	oldError = INT_error(auxiMat, n, n);

	cout << '\t' << oldError;
	
	double oldMax = maxValue(auxMat, n, n);

	oldImgQ = imgQ = imageQuality(auxMat, inputMat, oldMax, n);

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyMatrix(inputMat, auxMat, n);
		times[i] = test13_DevelopedDecompComp(auxMat, auxiMat, n, standard, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean) << '\t';
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyMatrix(inputiMat, auxiMat, n);
	test13_DevelopedDecompComp(auxMat, auxiMat, n, standard, false);

	error = INT_error(auxiMat, n, n);

	cout << '\t' << error;
	
	double max = maxValue(auxMat, n, n);

	imgQ = imageQuality(auxMat, inputMat, max, n);

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	cout << endl;
	cout << "Quality gain (%)" << endl;
	cout << endl;

	cout << "EUC" << "\t\t" << relativeGain(imgQ.euc, oldImgQ.euc) << endl;
	cout << "MSE" << "\t\t" << relativeGain(imgQ.mse, oldImgQ.mse) << endl;
	cout << "PSNR" << "\t\t" << -relativeGain(imgQ.psnr, oldImgQ.psnr) << endl;
	cout << "Error" << "\t\t" << relativeGain(error, oldError) << endl;

	cout << endl;

	delete [] times;
}

int test13(int argc, char **argv)
{
	if (argc < 1)
	{
		cout << endl;
		cout << "	ERROR: lack of parameters." << endl;
		test13Desc();
		test13Param();
		return 1;
	}

	uint n = atoi(argv[0]);

	if (n <= 0)
	{
		cout << endl;
		cout << "	ERROR: <order> has to be greater than 0." << endl;
		test13Desc();
		test13Param();
		return 1;
	}
	else if (!isPowerOfTwo(n))
	{
		cout << endl;
		cout << "	ERROR: <order> has to be power of two." << endl;
		test13Desc();
		test13Param();
		return 1;
	}

	test13Desc();
	cout << endl;
	cout << "Computing matrix of order " << n  << "..." << endl;
	cout << n << " X " << n << " = " << n * n << "..." << endl;
	cout << endl;

	double **inputMat = new double*[n];		// Matrices used as input for all executions.
	interval **inputiMat = new interval*[n];

	double **auxMat = new double*[n];		// Auxiliary matrices used for each execution.
	interval **auxiMat = new interval*[n];

	// Allocate and fill input matrices with n*n random values.
	for (uint i = 0; i < n; ++i)
	{
		inputMat[i] = new double[n];
		inputiMat[i] = new interval[n];

		auxMat[i] = new double[n];
		auxiMat[i] = new interval[n];

		for (uint j = 0; j < n; ++j)
		{
			inputMat[i][j] = rand() % n;
			inputiMat[i][j] = interval(inputMat[i][j]);
		}
	}

	test13_Script(inputMat, auxMat, inputiMat, auxiMat, n, true);
	test13_Script(inputMat, auxMat, inputiMat, auxiMat, n, false);

	// Deallocating memory.

	for (uint i = 0; i < n; ++i)
	{
		delete [] inputMat[i];
		delete [] inputiMat[i];
		delete [] auxMat[i];
		delete [] auxiMat[i];
	}

	delete [] inputMat;
	delete [] inputiMat;
	delete [] auxMat;
	delete [] auxiMat;

	return 0;
}

void test12Desc()
{
	cout << endl;
	cout << "	==================== Test 1 ====================" << endl;
	cout << endl;
	cout << "	* This test is about mesurement of time and calculus " << endl;
	cout << "	exactitude for the one-dimensional Haar Wavelet Transform." << endl;
}

void test12Param()
{
	cout << endl;
	cout << "	Parameters: <size>" << endl;
	cout << "	" << endl;
	cout << "	<size> must be an unsigned integer greater" << endl;
	cout << "	than 0 and also be power of two." << endl;
	cout << endl;
}

double test12_OriginalDecomp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_Decomposition(vec, n, true);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_Decomposition(ivec, n, true);
	}	

	return timer;
}

double test12_DevelopedDecomp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_Decomposition(vec, n, false);
		VinisNormalization(vec, n);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_Decomposition(ivec, n, false);
		INT_VinisNormalization(ivec, n);
	}	

	return timer;
}

double test12_OriginalComp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_Composition(vec, n, true);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_Composition(ivec, n, true);
	}	

	return timer;
}

double test12_DevelopedComp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		VinisNormalization(vec, n, true);
		Haar_Composition(vec, n, false);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_VinisNormalization(ivec, n, true);
		INT_Haar_Composition(ivec, n, false);
	}	

	return timer;
}

double test12_OriginalDecompComp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_Decomposition(vec, n, true);
		Haar_Composition(vec, n, true);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_Decomposition(ivec, n, true);
		INT_Haar_Composition(ivec, n, true);
	}	

	return timer;
}

double test12_DevelopedDecompComp(double *vec, interval *ivec, int n, bool getTime)
{
	double timer = 0.0;

	if (getTime)
	{
		// Get time of execution.
		startTimeCounter();
		Haar_Decomposition(vec, n, false);
		VinisNormalization(vec, n);
		VinisNormalization(vec, n, true);
		Haar_Composition(vec, n, false);
		timer = getTimeCounter();
	}
	else
	{
		// Interval transformation for error calculation.
		INT_Haar_Decomposition(ivec, n, false);
		INT_VinisNormalization(ivec, n);
		INT_VinisNormalization(ivec, n, true);
		INT_Haar_Composition(ivec, n, false);
	}	

	return timer;
}

int test12(int argc, char **argv)
{
	if (argc < 1)
	{
		cout << endl;
		cout << "	ERROR: lack of parameters." << endl;
		test12Desc();
		test12Param();
		return 1;
	}

	uint n = atoi(argv[0]);

	if (n <= 0)
	{
		cout << endl;
		cout << "	ERROR: <size> has to be greater than 0." << endl;
		test12Desc();
		test12Param();
		return 1;
	}
	else if (!isPowerOfTwo(n))
	{
		cout << endl;
		cout << "	ERROR: <size> has to be power of two." << endl;
		test12Desc();
		test12Param();
		return 1;
	}

	test12Desc();
	cout << endl;
	cout << "Computing vector of size " << n  << "..." << endl;
	cout << endl;

	real error, oldError;

	TimeMesurement timeMesure, oldMesure;

	double *inputVec = new double[n];		// Vectors used as input for all executions.
	interval *inputiVec = new interval[n];

	double *auxVec = new double[n];			// Auxiliary vectors used for each execution.
	interval *auxiVec = new interval[n];

	double *times = new double[30];			// Vector to store the time of each execution.

	// Fill input vectors with n random values.
	for (uint i = 0; i < n; ++i)
	{
		inputVec[i] = rand() %  n;
		inputiVec[i] = interval(inputVec[i]);
	}

	cout << "========== Decomposition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_OriginalDecomp(auxVec, auxiVec, n, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_OriginalDecomp(auxVec, auxiVec, n, false);

	oldError = INT_error(auxiVec, n);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_DevelopedDecomp(auxVec, auxiVec, n, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean);
	cout << "\t\t" << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_DevelopedDecomp(auxVec, auxiVec, n, false);

	error = INT_error(auxiVec, n);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	cout << endl;
	cout << "========== Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "Error (%)" << endl;
	
	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_OriginalComp(auxVec, auxiVec, n, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_OriginalComp(auxVec, auxiVec, n, false);

	oldError = INT_error(auxiVec, n);

	cout << '\t' << oldError;
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_DevelopedComp(auxVec, auxiVec, n, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean);
	cout << "\t\t" << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_DevelopedComp(auxVec, auxiVec, n, false);

	error = INT_error(auxiVec, n);

	cout << '\t' << error;
	cout << '\t' << relativeGain(error, oldError);
	cout << endl;

	ImageQuality<double> imgQ, oldImgQ;

	cout << endl;
	cout << "========== Decomposition & Composition: " << endl;
	cout << "\t\t";
	cout << "Speedup" << "\t\t";
	cout << "Time" << "\t\t";
	cout << "StdDev" << "\t\t";
	cout << "Error" << "\t\t";
	cout << "EUC" << "\t\t";
	cout << "MSE" << "\t\t";
	cout << "PSNR" << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_OriginalDecompComp(auxVec, auxiVec, n, true);
	}

	oldMesure = timeMesure = runTimeMesurement(times, 30);

	cout << "Original";
	cout << "\t\t";
	cout << '\t' << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_OriginalDecompComp(auxVec, auxiVec, n, false);

	oldError = INT_error(auxiVec, n);
	
	cout << '\t' << oldError;

	double oldMax = maxValue(auxVec, n);

	oldImgQ = imgQ = imageQuality(auxVec, inputVec, oldMax, n);

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	for (int i = 0; i < 30; ++i)
	{
		copyVector(inputVec, auxVec, n);
		times[i] = test12_DevelopedDecompComp(auxVec, auxiVec, n, true);
	}

	timeMesure = runTimeMesurement(times, 30);

	cout << "Developed";
	cout << '\t' << speedupCalc(timeMesure.mean, oldMesure.mean);
	cout << "\t\t" << timeMesure.mean;
	cout << '\t' << timeMesure.stdDev;

	copyVector(inputiVec, auxiVec, n);
	test12_DevelopedDecompComp(auxVec, auxiVec, n, false);

	error = INT_error(auxiVec, n);

	cout << '\t' << error;
	
	double max = maxValue(auxVec, n);

	imgQ = imageQuality(auxVec, inputVec, n, n);

	cout << '\t' << imgQ.euc;
	cout << '\t' << imgQ.mse;
	cout << '\t' << imgQ.psnr << endl;

	cout << endl;
	cout << "Quality gain (%)" << endl;
	cout << endl;

	cout << "EUC" << "\t\t" << relativeGain(imgQ.euc, oldImgQ.euc) << endl;
	cout << "MSE" << "\t\t" << relativeGain(imgQ.mse, oldImgQ.mse) << endl;
	cout << "PSNR" << "\t\t" << -relativeGain(imgQ.psnr, oldImgQ.psnr) << endl;
	cout << "Error" << "\t\t" << relativeGain(error, oldError) << endl;

	cout << endl;

	// Deallocating memory.
	delete [] inputVec;
	delete [] inputiVec;
	delete [] auxVec;
	delete [] auxiVec;
	delete [] times;

	return 0;

}

int test11(const char *filepath, int type)
{
	double **input = NULL;
	double **data = NULL;
	ImageQuality<double> imgQuality;
	EucMSE_data<double> imgMseEuc;
	interval **intInput = NULL;
	interval **intData = NULL;
	ImageQuality<interval> intImgQuality;
	EucMSE_data<interval> intImgMseEuc;
	ImageInfo imgInfo;

	input = carregar_imagem((char*)(filepath), &imgInfo);
	data = carregar_imagem((char*)(filepath), &imgInfo);

	if (input == NULL)
	{
		cout << "Error! Can't open file!" << endl;
		return 1;
	}

	intInput = new interval*[imgInfo.x];
	intData = new interval*[imgInfo.x];

	for (int i = 0; i < imgInfo.x; ++i)
	{
		intInput[i] = new interval[imgInfo.y];
		intData[i] = new interval[imgInfo.y];

		for (int j = 0; j < imgInfo.y; ++j)
		{
			intInput[i][j] = interval(input[i][j]);
			intData[i][j] = interval(data[i][j]);
		}
	}

	if (type == 0)
	{
		cout << "========== Original Standard Transform" << endl;

		Haar_MatrixDecomposition(data, imgInfo.x, imgInfo.y, true, true);
		Haar_MatrixComposition(data, imgInfo.x, imgInfo.y, true, true);

		INT_Haar_MatrixDecomposition(intData, imgInfo.x, imgInfo.y, true, true);
		INT_Haar_MatrixComposition(intData, imgInfo.x, imgInfo.y, true, true);
	}
	else if (type == 1)
	{
		cout << "========== New Standard Transform" << endl;

		Haar_MatrixDecomposition(data, imgInfo.x, imgInfo.y, false, true);
		VinisMatrixNormalization(data, imgInfo.x, true);
		VinisMatrixNormalization(data, imgInfo.x, true, true);
		Haar_MatrixComposition(data, imgInfo.x, imgInfo.y, false, true);

		INT_Haar_MatrixDecomposition(intData, imgInfo.x, imgInfo.y, false, true);
		INT_VinisMatrixNormalization(intData, imgInfo.x, true);
		INT_VinisMatrixNormalization(intData, imgInfo.x, true, true);
		INT_Haar_MatrixComposition(intData, imgInfo.x, imgInfo.y, false, true);
	}
	else if (type == 2)
	{
		cout << "========== Original Non Standard Transform" << endl;

		Haar_MatrixDecomposition(data, imgInfo.x, imgInfo.y, true, false);
		Haar_MatrixComposition(data, imgInfo.x, imgInfo.y, true, false);

		INT_Haar_MatrixDecomposition(intData, imgInfo.x, imgInfo.y, true, false);
		INT_Haar_MatrixComposition(intData, imgInfo.x, imgInfo.y, true, false);
	}
	else if (type == 3)
	{
		cout << "========== New Non Standard Transform" << endl;

		Haar_MatrixDecomposition(data, imgInfo.x, imgInfo.y, false, false);
		VinisMatrixNormalization(data, imgInfo.x, false);
		VinisMatrixNormalization(data, imgInfo.x, false, true);
		Haar_MatrixComposition(data, imgInfo.x, imgInfo.y, false, false);

		INT_Haar_MatrixDecomposition(intData, imgInfo.x, imgInfo.y, false, false);
		INT_VinisMatrixNormalization(intData, imgInfo.x, false);
		INT_VinisMatrixNormalization(intData, imgInfo.x, false, true);
		INT_Haar_MatrixComposition(intData, imgInfo.x, imgInfo.y, false, false);
	}

	imgMseEuc = EucMSE(data, input, imgInfo.x);
	imgQuality.euc = imgMseEuc.euc;
	imgQuality.mse = imgMseEuc.mse;
	imgQuality.psnr = PSNR(imgQuality.mse, 255.0);

	intImgMseEuc = INT_EucMSE(intData, intInput, imgInfo.x);
	intImgQuality.euc = intImgMseEuc.euc;
	intImgQuality.mse = intImgMseEuc.mse;
	intImgQuality.psnr = INT_PSNR(intImgQuality.mse, interval(255.0));

	//escrever_imagem("out.ppm", result, imgInfo);

	cout << "Error:" << '\t' << INT_error(intData, imgInfo.x, imgInfo.x) << endl;
	cout << "EUC:" << '\t' << imgMseEuc.euc << '\t' << intImgMseEuc.euc << endl;
	cout << "MSE: " << '\t' << imgMseEuc.mse << '\t' << intImgMseEuc.mse << endl;
	cout << "PSNR: " << '\t' << imgQuality.psnr << '\t' << '\t' << intImgQuality.psnr << endl;
	cout << endl;

	// Dealocates memory.
	for (int i = 0; i < imgInfo.x; ++i)
	{
		delete [] input[i];
		delete [] data[i];
		delete [] intInput[i];
		delete [] intData[i];
	}

	delete [] input;
	delete [] data;
	delete [] intInput;
	delete [] intData;

	return 0;
}

int test10(const char *filepath, int type, int levels)
{
	double **input = NULL;
	double ***data = NULL;
	double **result = NULL;
	ImageQuality<double> imgQuality;
	EucMSE_data<double> imgMseEuc;
	interval **intInput = NULL;
	interval ***intData = NULL;
	interval **intResult = NULL;
	ImageQuality<interval> intImgQuality;
	EucMSE_data<interval> intImgMseEuc;
	ImageInfo imgInfo;

	input = carregar_imagem((char*)(filepath), &imgInfo);

	if (input == NULL)
	{
		cout << "Error! Can't open file!" << endl;
		return 1;
	}

	intInput = new interval*[imgInfo.x];

	for (int i = 0; i < imgInfo.x; ++i)
	{
		intInput[i] = new interval[imgInfo.y];

		for (int j = 0; j < imgInfo.y; ++j)
		{
			intInput[i][j] = interval(input[i][j]);
		}
	}

	if (type == 0)
	{
		cout << "========== Original Standard Transform" << endl;

		data = Haar_atrous_MatrixDecomposition(input, imgInfo.x, imgInfo.y, levels, true, true);
		result = Haar_atrous_MatrixComposition(data, imgInfo.x, imgInfo.y, levels, true);

		intData = Haar_atrous_MatrixDecomposition(intInput, imgInfo.x, imgInfo.y, levels, true, true);
		intResult = Haar_atrous_MatrixComposition(intData, imgInfo.x, imgInfo.y, levels, true);
	}
	else if (type == 1)
	{
		cout << "========== New Standard Transform" << endl;

		data = Haar_atrous_MatrixDecomposition(input, imgInfo.x, imgInfo.y, levels, false, true);
		Haar_atrous_MatrixNormalization(input, data, imgInfo.x, imgInfo.y, levels);
		result = Haar_atrous_MatrixComposition(data, imgInfo.x, imgInfo.y, levels, true);

		intData = Haar_atrous_MatrixDecomposition(intInput, imgInfo.x, imgInfo.y, levels, false, true);
		Haar_atrous_MatrixNormalization(intInput, intData, imgInfo.x, imgInfo.y, levels);
		intResult = Haar_atrous_MatrixComposition(intData, imgInfo.x, imgInfo.y, levels, true);
	}
	else if (type == 2)
	{
		cout << "========== Original Non Standard Transform" << endl;

		data = Haar_atrous_MatrixDecomposition(input, imgInfo.x, imgInfo.y, levels, true, false);
		result = Haar_atrous_MatrixComposition(data, imgInfo.x, imgInfo.y, levels, false);

		intData = Haar_atrous_MatrixDecomposition(intInput, imgInfo.x, imgInfo.y, levels, true, false);
		intResult = Haar_atrous_MatrixComposition(intData, imgInfo.x, imgInfo.y, levels, false);
	}
	else if (type == 3)
	{
		cout << "========== New Non Standard Transform" << endl;

		data = Haar_atrous_MatrixDecomposition(input, imgInfo.x, imgInfo.y, levels, false, false);
		Haar_atrous_MatrixNormalization(input, data, imgInfo.x, imgInfo.y, levels);
		result = Haar_atrous_MatrixComposition(data, imgInfo.x, imgInfo.y, levels, false);

		intData = Haar_atrous_MatrixDecomposition(intInput, imgInfo.x, imgInfo.y, levels, false, false);
		Haar_atrous_MatrixNormalization(intInput, intData, imgInfo.x, imgInfo.y, levels);
		intResult = Haar_atrous_MatrixComposition(intData, imgInfo.x, imgInfo.y, levels, false);
	}

	imgMseEuc = EucMSE(result, input, imgInfo.x);
	imgQuality.euc = imgMseEuc.euc;
	imgQuality.mse = imgMseEuc.mse;
	imgQuality.psnr = PSNR(imgQuality.mse, 255.0);

	intImgMseEuc = INT_EucMSE(intResult, intInput, imgInfo.x);
	intImgQuality.euc = intImgMseEuc.euc;
	intImgQuality.mse = intImgMseEuc.mse;
	intImgQuality.psnr = INT_PSNR(intImgQuality.mse, interval(255.0));

	//escrever_imagem("out.ppm", result, imgInfo);

	cout << "Error:" << '\t' << INT_error(intResult, imgInfo.x, imgInfo.x) << endl;
	cout << "EUC:" << '\t' << imgMseEuc.euc << '\t' << intImgMseEuc.euc << endl;
	cout << "MSE: " << '\t' << imgMseEuc.mse << '\t' << intImgMseEuc.mse << endl;
	cout << "PSNR: " << '\t' << imgQuality.psnr << '\t' << '\t' << intImgQuality.psnr << endl;
	cout << endl;

	// Dealocates memory.
	for (int i = 0; i < imgInfo.x; ++i)
	{
		delete [] input[i];
		delete [] intInput[i];
	}

	delete [] input;
	delete [] intInput;

	return 0;
}

int test9(int n, int max, int levels)
{
	interval **intMatrix = new interval*[n];
	interval **intMat = new interval*[n];
	interval ***intData;
	ImageQuality<interval> imgQuality;
	EucMSE_data<interval> imgMseEuc;

	for (int i = 0; i < n; i++)
	{
		intMatrix[i] = new interval[n];
		intMat[i] = new interval[n];

		for (int j = 0; j < n; j++)
		{
			intMatrix[i][j] = interval(((i + 1) + (j + 1)) % max);
			intMat[i][j] = interval(((i + 1) + (j + 1)) % max);
		}
	}

	cout << "==================== Original Algorithms" << endl;

	copyMatrix(intMatrix, intMat, n);
	intData = Haar_atrous_MatrixDecomposition(intMat, n, n, levels, true, true);
	deleteMatrix<interval>(intMat, n);
	intMat = Haar_atrous_MatrixComposition(intData, n, n, levels, true);

	//imgQuality = imageQuality<interval>(intMat, intMatrix, max, n);
	imgMseEuc = INT_EucMSE(intMat, intMatrix, n);
	imgQuality.euc = imgMseEuc.euc;
	imgQuality.mse = imgMseEuc.mse;
	imgQuality.psnr = INT_PSNR(imgQuality.mse, interval(max));

	cout << "---------- Standard Normalized:" << '\t' << INT_error(intMat, n, n) << endl;
	cout << "---------- EUC: " << imgMseEuc.euc << endl;
	cout << "---------- MSE: " << imgMseEuc.mse << endl;
	cout << "---------- PSNR: " << imgQuality.psnr << endl;
	cout << endl;


	return 0;
}

int test8(int n, int levels)
{
	if (n <= 0 || levels <= 0) return 1;

	interval *intVector = new interval[n];
	interval *intVec = new interval[n];
	interval **intData;

	for (int i = 0; i < n; ++i)
	{
		intVector[i] = interval((rand() %  2048) - 1024);
	}

	cout << "==================== Original Algorithms" << endl;

	copyVector(intVector, intVec, n);
	intData = Haar_atrous_Decomposition(intVec, n, levels, false);
	delete [] intVec;
	intVec = Haar_atrous_Composition(intData, n, levels);

	cout << "---------- Non-Normalized:" << '\t' << INT_error(intVec, n) << endl;

	copyVector(intVector, intVec, n);
	intData = Haar_atrous_Decomposition(intVec, n, levels, true);
	delete [] intVec;
	intVec = Haar_atrous_Composition(intData, n, levels);

	cout << "---------- Normalized:" << "\t\t" << INT_error(intVec, n) << endl;

	cout << "==================== Our New Algorithms" << endl;

	copyVector(intVector, intVec, n);
	intData = Haar_atrous_Decomposition(intVec, n, levels, false);
	delete [] intVec;
	intVec = Haar_atrous_Composition(intData, n, levels);

	cout << "---------- Non-Normalized:" << '\t' << INT_error(intVec, n) << endl;

	copyVector(intVector, intVec, n);
	intData = Haar_atrous_Decomposition(intVec, n, levels, false);
	Haar_atrous_Normalization(intVec, intData, n, levels);
	Haar_atrous_Normalization(intVec, intData, n, levels, true);
	delete [] intVec;
	intVec = Haar_atrous_Composition(intData, n, levels);

	cout << "---------- Normalized:" << "\t\t" << INT_error(intVec, n) << endl;

	return 0;
}

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
	cout << "========= Using original algorithms:" << endl;
	cout << endl;
	cout << "* vector states by row!" << endl;
	cout << endl;

	printVectors<double>(resultBuffer, 8, 4, 10);
	cout << endl;
	cout << "* vector states by row!" << endl;
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
	cout << "========= Using new algorithms:" << endl;
	cout << endl;
	cout << "* vector states by row!" << endl;
	cout << endl;

	printVectors<double>(resultBuffer, 8, 4, 10);
	cout << endl;
	cout << "* vector states by row!" << endl;
	cout << endl;
	printVectors<interval>(intResultBuffer, 8, 4, 1);
	
	// =====================================================
	// ==================== Matrix part ====================
	// =====================================================

	double **mat = new double*[8];
	double **matrix = new double*[8];
	double ***resultMatBuffer = new double**[4];
	interval **intMat = new interval*[8];
	interval **intMatrix = new interval*[8];
	interval ***intResultMatBuffer = new interval**[4];

	for (int i = 0; i < 8; ++i)
	{
		mat[i] = new double[8];
		matrix[i] = new double[8];
		intMat[i] = new interval[8];
		intMatrix[i] = new interval[8];
	}
	
	for (int i = 0; i < 8; ++i)
	for (int j = 0; j < 8; ++j)
	{
		matrix[i][j] = (double)(i + j);
		intMatrix[i][j] = interval(i + j);
	}
	
	for (int i = 0; i < 4; ++i)
	{
		resultMatBuffer[i] = new double*[8];
		intResultMatBuffer[i] = new interval*[8];

		for (int j = 0; j < 8; ++j)
		{
			resultMatBuffer[i][j] = new double[8];
			intResultMatBuffer[i][j] = new interval[8];
		}
	}
	
	copyMatrix<double>(matrix, mat, 8);
	copyMatrix<interval>(intMatrix, intMat, 8);
	
	copyMatrix<double>(mat, resultMatBuffer[0], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[0], 8);
	
	Haar_MatrixDecomposition(mat, 8, 8, true, true);
	INT_Haar_MatrixDecomposition(intMat, 8, 8, true, true);
	
	copyMatrix<double>(mat, resultMatBuffer[1], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[1], 8);
	
	Haar_Matrix_Compression(mat, 8, 0.01);
	INT_Haar_Matrix_Compression(intMat, 8, 0.01);
	
	copyMatrix<double>(mat, resultMatBuffer[2], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[2], 8);
	
	Haar_MatrixComposition(mat, 8, 8, true, true);
	INT_Haar_MatrixComposition(intMat, 8, 8, true, true);
	
	copyMatrix<double>(mat, resultMatBuffer[3], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[3], 8);
	
	cout << endl;
	cout << "=====================================================" << endl;
	cout << "==================== Matrix part ====================" << endl;
	cout << "=====================================================" << endl;
	cout << endl;
	cout << "========== Using Original Standard Transformation" << endl;
	cout << endl;
	
	printMatrices<double>(resultMatBuffer, 8, 8, 4, 10);
	cout << endl;
	printMatrices<interval>(intResultMatBuffer, 8, 8, 4, 1);
	
	copyMatrix<double>(matrix, mat, 8);
	copyMatrix<interval>(intMatrix, intMat, 8);
	
	copyMatrix<double>(mat, resultMatBuffer[0], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[0], 8);

	Haar_MatrixDecomposition(mat, 8, 8, true, false);
	INT_Haar_MatrixDecomposition(intMat, 8, 8, true, false);
	
	copyMatrix<double>(mat, resultMatBuffer[1], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[1], 8);
	
	Haar_Matrix_Compression(mat, 8, 0.01);
	INT_Haar_Matrix_Compression(intMat, 8, 0.01);
	
	copyMatrix<double>(mat, resultMatBuffer[2], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[2], 8);
	
	Haar_MatrixComposition(mat, 8, 8, true, false);
	INT_Haar_MatrixComposition(intMat, 8, 8, true, false);
	
	copyMatrix<double>(mat, resultMatBuffer[3], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[3], 8);

	cout << endl;
	cout << "========== Using Original Non-Standard Transformation" << endl;
	cout << endl;

	printMatrices<double>(resultMatBuffer, 8, 8, 4, 10);
	cout << endl;
	printMatrices<interval>(intResultMatBuffer, 8, 8, 4, 1);

	copyMatrix<double>(matrix, mat, 8);
	copyMatrix<interval>(intMatrix, intMat, 8);
	
	copyMatrix<double>(mat, resultMatBuffer[0], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[0], 8);

	Haar_MatrixDecomposition(mat, 8, 8, false, true);
	VinisMatrixNormalization(mat, 8, true);
	INT_Haar_MatrixDecomposition(intMat, 8, 8, false, true);
	INT_VinisMatrixNormalization(intMat, 8, true);
	
	copyMatrix<double>(mat, resultMatBuffer[1], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[1], 8);
	
	Haar_Matrix_Compression(mat, 8, 0.01);
	INT_Haar_Matrix_Compression(intMat, 8, 0.01);
	
	copyMatrix<double>(mat, resultMatBuffer[2], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[2], 8);
	
	VinisMatrixNormalization(mat, 8, true, true);
	Haar_MatrixComposition(mat, 8, 8, false, true);
	INT_VinisMatrixNormalization(intMat, 8, true, true);
	INT_Haar_MatrixComposition(intMat, 8, 8, false, true);
	
	copyMatrix<double>(mat, resultMatBuffer[3], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[3], 8);

	cout << endl;
	cout << "========== Using Our New Standard Transformation" << endl;
	cout << endl;

	printMatrices<double>(resultMatBuffer, 8, 8, 4, 10);
	cout << endl;
	printMatrices<interval>(intResultMatBuffer, 8, 8, 4, 1);

	copyMatrix<double>(matrix, mat, 8);
	copyMatrix<interval>(intMatrix, intMat, 8);
	
	copyMatrix<double>(mat, resultMatBuffer[0], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[0], 8);

	Haar_MatrixDecomposition(mat, 8, 8, false, false);
	VinisMatrixNormalization(mat, 8, false);
	INT_Haar_MatrixDecomposition(intMat, 8, 8, false, false);
	INT_VinisMatrixNormalization(intMat, 8, false);
	
	copyMatrix<double>(mat, resultMatBuffer[1], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[1], 8);
	
	Haar_Matrix_Compression(mat, 8, 0.01);
	INT_Haar_Matrix_Compression(intMat, 8, 0.01);
	
	copyMatrix<double>(mat, resultMatBuffer[2], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[2], 8);
	
	VinisMatrixNormalization(mat, 8, false, true);
	Haar_MatrixComposition(mat, 8, 8, false, false);
	INT_VinisMatrixNormalization(intMat, 8, false, true);
	INT_Haar_MatrixComposition(intMat, 8, 8, false, false);
	
	copyMatrix<double>(mat, resultMatBuffer[3], 8);
	copyMatrix<interval>(intMat, intResultMatBuffer[3], 8);

	cout << endl;
	cout << "========== Using Our New Non-Standard Transformation" << endl;
	cout << endl;

	printMatrices<double>(resultMatBuffer, 8, 8, 4, 10);
	cout << endl;
	printMatrices<interval>(intResultMatBuffer, 8, 8, 4, 1);

}

void fundamentalTest2()
{
	double *vec = new double[8];
	double **data, *result;
	interval *intVec = new interval[8];
	interval **intData, *intResult;

	vec[0] = 1.0; vec[1] = 7.0; vec[2] = 2.0; vec[3] = 6.0;
	vec[4] = 3.0; vec[5] = 5.0; vec[6] = 4.0; vec[7] = 4.0;

	for (int i = 0; i < 8; ++i) intVec[i] = interval(vec[i]);

	data = Haar_atrous_Decomposition<double>(vec, 8, 4, true);
	result =  Haar_atrous_Composition<double>(data, 8, 4);
	intData = Haar_atrous_Decomposition<interval>(intVec, 8, 4, true);
	intResult = Haar_atrous_Composition<interval>(intData, 8, 4);
	
	cout << endl;
	cout << "========== Normalized Decomposition" << endl;
	cout << endl;
	cout << endl;

	cout << "---------- Input" << endl;
	printVectors<double>(&vec, 8, 1, 10);
	cout << endl;
	cout << "---------- Normalized Decomposition" << endl;
	printVectors<double>(data, 8, 8, 10);
	cout << endl;
	cout << "---------- Normalized Composition" << endl;
	printVectors<double>(&result, 8, 1, 10);
	cout << endl;
	cout << "---------- Input" << endl;
	printVectors<interval>(&intVec, 8, 1, 1);
	cout << endl;
	cout << "---------- Normalized Decomposition" << endl;
	printVectors<interval>(intData, 8, 8, 1);
	cout << "---------- Normalized Composition" << endl;
	printVectors<interval>(&intResult, 8, 1, 1);
	cout << endl;

	cout << endl;
	cout << "========== Normalized Decomposition using our new methods" << endl;
	cout << endl;
	cout << endl;

	deleteMatrix<double>(data, 8);
	deleteMatrix<interval>(intData, 8);

	data = Haar_atrous_Decomposition<double>(vec, 8, 4, false);

	cout << "---------- Input" << endl;
	printVectors<double>(&vec, 8, 1, 10);
	cout << endl;
	cout << "---------- Non-Normalized Decomposition" << endl;
	printVectors<double>(data, 8, 8, 10);
	cout << endl;

	Haar_atrous_Normalization<double>(vec, data, 8, 4);

	cout << "---------- Normalization Step" << endl;
	printVectors<double>(data, 8, 8, 10);
	cout << endl;

	Haar_atrous_Normalization<double>(vec, data, 8, 4, true);

	cout << "---------- Denormalization Step" << endl;
	printVectors<double>(data, 8, 8, 10);
	cout << endl;

	delete [] result;

	result = Haar_atrous_Composition<double>(data, 8, 4);

	cout << "---------- Non-Normalized Composition" << endl;
	printVectors<double>(&result, 8, 1, 10);
	cout << endl;

	intData = Haar_atrous_Decomposition<interval>(intVec, 8, 4, false);

	cout << "---------- Input" << endl;
	printVectors<interval>(&intVec, 8, 1, 1);
	cout << endl;
	cout << "---------- Non-Normalized Decomposition" << endl;
	printVectors<interval>(intData, 8, 8, 1);
	cout << endl;

	Haar_atrous_Normalization<interval>(intVec, intData, 8, 4);

	cout << "---------- Normalization Step" << endl;
	printVectors<interval>(intData, 8, 8, 1);
	cout << endl;

	Haar_atrous_Normalization<interval>(intVec, intData, 8, 4, true);

	cout << "---------- Denormalization Step" << endl;
	printVectors<interval>(intData, 8, 8, 1);
	cout << endl;

	delete [] intResult;

	intResult = Haar_atrous_Composition<interval>(intData, 8, 4);

	cout << "---------- Non-Normalized Composition" << endl;
	printVectors<interval>(&intResult, 8, 1, 1);
	cout << endl;

}

void fundamentalTest3()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (4) test is about non-normalized Haar A-Trous " << endl;
	cout << "* decomposition and composition. It uses a small 4x4 matrix in order " << endl;
	cout << "* to help the visualization of the result and further comparison with " << endl;
	cout << "* its concept." << endl;
	cout << "********************" << endl;
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
	
	double **matrix = new double*[4];
	
	for (int m = 0; m < 4; m++)
		matrix[m] = new double[4];
	
	copyMatrix<double>(mat, matrix, 4);

	interval **intMatrix = new interval*[4];
	
	for (int m = 0; m < 4; m++)
		intMatrix[m] = new interval[4];
	
	copyMatrix<interval>(intmat, intMatrix, 4);

	double ***data = NULL;
	interval ***intData = NULL;

	double **restore;
	interval **intRestore;

	cout << endl;
	cout << "// ==================== Standard Transformation" << endl;
	cout << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Input: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, false, true);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, true);
	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, false, true);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, true);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

	for (int i = 0; i < 4 * 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			delete [] data[i][j];
			delete [] intData[i][j];
		}

		delete [] data[i];
		delete [] intData[i];
	}

	delete [] data;
	delete [] intData;

	copyMatrix<double>(mat, matrix, 4);
	copyMatrix<interval>(intmat, intMatrix, 4);

	cout << endl;
	cout << "// ==================== Non Standard Transformation" << endl;
	cout << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Input: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, false, false);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, false);
	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, false, false);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, false);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

}

void fundamentalTest4()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (4) test is about normalized Haar A-Trous " << endl;
	cout << "* decomposition and composition. It uses a small 4x4 matrix in order " << endl;
	cout << "* to help the visualization of the result and further comparison with " << endl;
	cout << "* its concept." << endl;
	cout << "********************" << endl;
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
	
	double **matrix = new double*[4];
	
	for (int m = 0; m < 4; m++)
		matrix[m] = new double[4];
	
	copyMatrix<double>(mat, matrix, 4);

	interval **intMatrix = new interval*[4];
	
	for (int m = 0; m < 4; m++)
		intMatrix[m] = new interval[4];
	
	copyMatrix<interval>(intmat, intMatrix, 4);

	double ***data = NULL;
	interval ***intData = NULL;

	double **restore;
	interval **intRestore;

	cout << endl;
	cout << "// ==================== Original Standard Transformation" << endl;
	cout << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Input: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, true, true);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, true);
	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, true, true);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, true);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

	for (int i = 0; i < 4 * 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			delete [] data[i][j];
			delete [] intData[i][j];
		}

		delete [] data[i];
		delete [] intData[i];
	}

	delete [] data;
	delete [] intData;

	cout << endl;
	cout << "// ==================== New Standard Transformation" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, false, true);
	Haar_atrous_MatrixNormalization(matrix, data, 4, 4, 2);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, true);

	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, false, true);
	Haar_atrous_MatrixNormalization(intMatrix, intData, 4, 4, 2);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, true);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

	for (int i = 0; i < 4 * 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			delete [] data[i][j];
			delete [] intData[i][j];
		}

		delete [] data[i];
		delete [] intData[i];
	}

	delete [] data;
	delete [] intData;

	copyMatrix<double>(mat, matrix, 4);
	copyMatrix<interval>(intmat, intMatrix, 4);

	cout << endl;
	cout << "// ==================== Original Non Standard Transformation" << endl;
	cout << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Input: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, true, false);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, false);
	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, true, false);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, false);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

	for (int i = 0; i < 4 * 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			delete [] data[i][j];
			delete [] intData[i][j];
		}

		delete [] data[i];
		delete [] intData[i];
	}

	delete [] data;
	delete [] intData;

	copyMatrix<double>(mat, matrix, 4);
	copyMatrix<interval>(intmat, intMatrix, 4);

	cout << endl;
	cout << "// ==================== New Non Standard Transformation" << endl;
	cout << endl;

	cout << endl;
	
	cout << "==========================================" << endl;
	cout << "Input: \n" << endl;
	cout << "{5, 6, 1, 2}" << endl;
	cout << "{4, 2, 5, 5}" << endl;
	cout << "{3, 1, 7, 1}" << endl;
	cout << "{6, 3, 5, 1}" << endl;
	cout << endl;

	data = Haar_atrous_MatrixDecomposition(matrix, 4, 4, 2, false, false);
	Haar_atrous_MatrixNormalization(matrix, data, 4, 4, 2);
	restore = Haar_atrous_MatrixComposition(data, 4, 4, 2, false);

	intData = Haar_atrous_MatrixDecomposition(intMatrix, 4, 4, 2, false, false);
	Haar_atrous_MatrixNormalization(intMatrix, intData, 4, 4, 2);
	intRestore = Haar_atrous_MatrixComposition(intData, 4, 4, 2, false);

	cout << "=============== Punctual Process:" << endl;
	cout << endl;

	printMatrices(data, 4, 4, 8, 10);

	cout << endl;
	cout << "========== Punctual Result of Composition:" << endl;
	cout << endl;

	printMatrix(restore, 4);

	cout << endl;
	cout << "=============== Interval Process:" << endl;
	cout << endl;

	printMatrices(intData, 4, 4, 8, 1);

	cout << endl;
	cout << "========== Interval Result of Composition:" << endl;
	cout << endl;

	printMatrix(intRestore, 4);

	cout << endl;

	for (int i = 0; i < 4 * 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			delete [] data[i][j];
			delete [] intData[i][j];
		}

		delete [] data[i];
		delete [] intData[i];
	}

	delete [] data;
	delete [] intData;

}

void fundamentalTest5()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (5) test is about the normalized Daubechies (db2)" << endl;
	cout << "* decomposition and composition steps. It uses a 16 long vector in order " << endl;
	cout << "* to help the visualization of the result and further comparison with " << endl;
	cout << "* its concept." << endl;
	cout << "********************" << endl;
	cout << endl;

	double *vec = new double[16];

	vec[0] = 9.0;
	vec[1] = 7.0;
	vec[2] = 3.0;
	vec[3] = 5.0;
	vec[4] = 11.0;
	vec[5] = 6.0;
	vec[6] = 2.0;
	vec[7] = 1.0;
	vec[8] = 9.0;
	vec[9] = 7.0;
	vec[10] = 4.0;
	vec[11] = 3.0;
	vec[12] = 1.0;
	vec[13] = 6.0;
	vec[14] = 7.0;
	vec[15] = 2.0;

	cout << "Input:" << endl;
	cout << endl;
	printVector(vec, 16);
	cout << endl;
	cout << "From the literature:" << endl;
	cout << endl;

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	Daub_Decomposition(vec, 16, true);
	printVector(vec, 16);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	Daub_Composition(vec, 16, true);
	printVector(vec, 16);

	cout << endl;
	cout << "Developed:" << endl;
	cout << endl;

	vec[0] = 9.0;
	vec[1] = 7.0;
	vec[2] = 3.0;
	vec[3] = 5.0;
	vec[4] = 11.0;
	vec[5] = 6.0;
	vec[6] = 2.0;
	vec[7] = 1.0;
	vec[8] = 9.0;
	vec[9] = 7.0;
	vec[10] = 4.0;
	vec[11] = 3.0;
	vec[12] = 1.0;
	vec[13] = 6.0;
	vec[14] = 7.0;
	vec[15] = 2.0;

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	Daub_Decomposition(vec, 16, false);
	printVector(vec, 16);

	cout << "Normalization Step:" << endl;
	cout << endl;
	Daub_Normalization(vec, 16);
	printVector(vec, 16);
	
	cout << "Normalized Composition:" << endl;
	cout << endl;
	Daub_Composition(vec, 16, true);
	printVector(vec, 16);
	

	cout << endl;

	delete [] vec;
}

void fundamentalTest6()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (6) test is about the standard normalized Daubechies" << endl;
	cout << "* (db2) decomposition and composition steps. It uses a 8x8 matrix in order " << endl;
	cout << "* to help the visualization of the result and further comparison with " << endl;
	cout << "* its concept." << endl;
	cout << "********************" << endl;
	cout << endl;

	double **mat = new double*[8];

	for (int m = 0; m < 8; m++)
		mat[m] = new double[8];

	mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2; mat[0][4] = 1; mat[0][5] = 3; mat[0][6] = 4; mat[0][7] = 7;
	mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5; mat[1][4] = 4; mat[1][5] = 2; mat[1][6] = 1; mat[1][7] = 5;
	mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1; mat[2][4] = 6; mat[2][5] = 4; mat[2][6] = 9; mat[2][7] = 0;
	mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1; mat[3][4] = 1; mat[3][5] = 5; mat[3][6] = 3; mat[3][7] = 7;
	mat[4][0] = 5; mat[4][1] = 6; mat[4][2] = 1; mat[4][3] = 2; mat[4][4] = 3; mat[4][5] = 6; mat[4][6] = 4; mat[4][7] = 5;
	mat[5][0] = 4; mat[5][1] = 2; mat[5][2] = 5; mat[5][3] = 5; mat[5][4] = 8; mat[5][5] = 3; mat[5][6] = 2; mat[5][7] = 1;
	mat[6][0] = 3; mat[6][1] = 1; mat[6][2] = 7; mat[6][3] = 1; mat[6][4] = 7; mat[6][5] = 9; mat[6][6] = 5; mat[6][7] = 8;
	mat[7][0] = 6; mat[7][1] = 3; mat[7][2] = 5; mat[7][3] = 1; mat[7][4] = 5; mat[7][5] = 6; mat[7][6] = 1; mat[7][7] = 2;

	cout << "Input:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	cout << "==================== Original method (literature):" << endl;
	cout << endl;

	Daub_StandardDecomposition(mat, 8, 8, true);

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_StandardComposition(mat, 8, 8, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2; mat[0][4] = 1; mat[0][5] = 3; mat[0][6] = 4; mat[0][7] = 7;
	mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5; mat[1][4] = 4; mat[1][5] = 2; mat[1][6] = 1; mat[1][7] = 5;
	mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1; mat[2][4] = 6; mat[2][5] = 4; mat[2][6] = 9; mat[2][7] = 0;
	mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1; mat[3][4] = 1; mat[3][5] = 5; mat[3][6] = 3; mat[3][7] = 7;
	mat[4][0] = 5; mat[4][1] = 6; mat[4][2] = 1; mat[4][3] = 2; mat[4][4] = 3; mat[4][5] = 6; mat[4][6] = 4; mat[4][7] = 5;
	mat[5][0] = 4; mat[5][1] = 2; mat[5][2] = 5; mat[5][3] = 5; mat[5][4] = 8; mat[5][5] = 3; mat[5][6] = 2; mat[5][7] = 1;
	mat[6][0] = 3; mat[6][1] = 1; mat[6][2] = 7; mat[6][3] = 1; mat[6][4] = 7; mat[6][5] = 9; mat[6][6] = 5; mat[6][7] = 8;
	mat[7][0] = 6; mat[7][1] = 3; mat[7][2] = 5; mat[7][3] = 1; mat[7][4] = 5; mat[7][5] = 6; mat[7][6] = 1; mat[7][7] = 2;

	cout << "==================== Developed method:" << endl;
	cout << endl;

	Daub_StandardDecomposition(mat, 8, 8, false);

	cout << "Non Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_StandardNormalization(mat, 8);

	cout << "Normalization Step Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_StandardComposition(mat, 8, 8, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

}

void fundamentalTest7()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (7) test is about the non standard normalized Daubechies" << endl;
	cout << "* (db2) decomposition and composition steps. It uses a 8x8 matrix in order " << endl;
	cout << "* to help the visualization of the result and further comparison with " << endl;
	cout << "* its concept." << endl;
	cout << "********************" << endl;
	cout << endl;

	double **mat = new double*[8];

	for (int m = 0; m < 8; m++)
		mat[m] = new double[8];

	mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2; mat[0][4] = 1; mat[0][5] = 3; mat[0][6] = 4; mat[0][7] = 7;
	mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5; mat[1][4] = 4; mat[1][5] = 2; mat[1][6] = 1; mat[1][7] = 5;
	mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1; mat[2][4] = 6; mat[2][5] = 4; mat[2][6] = 9; mat[2][7] = 0;
	mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1; mat[3][4] = 1; mat[3][5] = 5; mat[3][6] = 3; mat[3][7] = 7;
	mat[4][0] = 5; mat[4][1] = 6; mat[4][2] = 1; mat[4][3] = 2; mat[4][4] = 3; mat[4][5] = 6; mat[4][6] = 4; mat[4][7] = 5;
	mat[5][0] = 4; mat[5][1] = 2; mat[5][2] = 5; mat[5][3] = 5; mat[5][4] = 8; mat[5][5] = 3; mat[5][6] = 2; mat[5][7] = 1;
	mat[6][0] = 3; mat[6][1] = 1; mat[6][2] = 7; mat[6][3] = 1; mat[6][4] = 7; mat[6][5] = 9; mat[6][6] = 5; mat[6][7] = 8;
	mat[7][0] = 6; mat[7][1] = 3; mat[7][2] = 5; mat[7][3] = 1; mat[7][4] = 5; mat[7][5] = 6; mat[7][6] = 1; mat[7][7] = 2;

	//// ********************** Comparison standard and non standard
	//
	//Daub_StandardDecomposition(mat, 8, 8, true);
	//cout << endl;
	//printMatrix(mat, 8);
	//cout << endl;
	//
	//mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2; mat[0][4] = 1; mat[0][5] = 3; mat[0][6] = 4; mat[0][7] = 7;
	//mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5; mat[1][4] = 4; mat[1][5] = 2; mat[1][6] = 1; mat[1][7] = 5;
	//mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1; mat[2][4] = 6; mat[2][5] = 4; mat[2][6] = 9; mat[2][7] = 0;
	//mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1; mat[3][4] = 1; mat[3][5] = 5; mat[3][6] = 3; mat[3][7] = 7;
	//mat[4][0] = 5; mat[4][1] = 6; mat[4][2] = 1; mat[4][3] = 2; mat[4][4] = 3; mat[4][5] = 6; mat[4][6] = 4; mat[4][7] = 5;
	//mat[5][0] = 4; mat[5][1] = 2; mat[5][2] = 5; mat[5][3] = 5; mat[5][4] = 8; mat[5][5] = 3; mat[5][6] = 2; mat[5][7] = 1;
	//mat[6][0] = 3; mat[6][1] = 1; mat[6][2] = 7; mat[6][3] = 1; mat[6][4] = 7; mat[6][5] = 9; mat[6][6] = 5; mat[6][7] = 8;
	//mat[7][0] = 6; mat[7][1] = 3; mat[7][2] = 5; mat[7][3] = 1; mat[7][4] = 5; mat[7][5] = 6; mat[7][6] = 1; mat[7][7] = 2;
	//
	//Daub_NonStandardDecomposition(mat, 8, 8, true);
	//cout << endl;
	//printMatrix(mat, 8);
	//cout << endl;
	//
	//return;
	//
	//// **********************

	cout << "Input:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	cout << "==================== Original method (literature):" << endl;
	cout << endl;

	Daub_NonStandardDecomposition(mat, 8, 8, true);

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_NonStandardComposition(mat, 8, 8, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	mat[0][0] = 5; mat[0][1] = 6; mat[0][2] = 1; mat[0][3] = 2; mat[0][4] = 1; mat[0][5] = 3; mat[0][6] = 4; mat[0][7] = 7;
	mat[1][0] = 4; mat[1][1] = 2; mat[1][2] = 5; mat[1][3] = 5; mat[1][4] = 4; mat[1][5] = 2; mat[1][6] = 1; mat[1][7] = 5;
	mat[2][0] = 3; mat[2][1] = 1; mat[2][2] = 7; mat[2][3] = 1; mat[2][4] = 6; mat[2][5] = 4; mat[2][6] = 9; mat[2][7] = 0;
	mat[3][0] = 6; mat[3][1] = 3; mat[3][2] = 5; mat[3][3] = 1; mat[3][4] = 1; mat[3][5] = 5; mat[3][6] = 3; mat[3][7] = 7;
	mat[4][0] = 5; mat[4][1] = 6; mat[4][2] = 1; mat[4][3] = 2; mat[4][4] = 3; mat[4][5] = 6; mat[4][6] = 4; mat[4][7] = 5;
	mat[5][0] = 4; mat[5][1] = 2; mat[5][2] = 5; mat[5][3] = 5; mat[5][4] = 8; mat[5][5] = 3; mat[5][6] = 2; mat[5][7] = 1;
	mat[6][0] = 3; mat[6][1] = 1; mat[6][2] = 7; mat[6][3] = 1; mat[6][4] = 7; mat[6][5] = 9; mat[6][6] = 5; mat[6][7] = 8;
	mat[7][0] = 6; mat[7][1] = 3; mat[7][2] = 5; mat[7][3] = 1; mat[7][4] = 5; mat[7][5] = 6; mat[7][6] = 1; mat[7][7] = 2;

	cout << "==================== Developed method:" << endl;
	cout << endl;

	Daub_NonStandardDecomposition(mat, 8, 8, false);

	cout << "Non Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_NonStandardNormalization(mat, 8);

	cout << "Normalization Step Decomposition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

	Daub_NonStandardComposition(mat, 8, 8, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat, 8);
	cout << endl;

}

void fundamentalTest8()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (8) test is about the standard normalized Daubechies" << endl;
	cout << "* (db2) decomposition and composition procedures. It uses a 16x16 matrix and" << endl;
	cout << "* performs 2 normalization steps for each procedure." << endl;
	cout << "********************" << endl;
	cout << endl;

	uint n = 16;

	double **mat1 = new double*[n];
	double **mat2 = new double*[n];

	for (uint i = 0; i < n; i++)
	{
		mat1[i] = new double[n];
		mat2[i] = new double[n];

		for (uint j = 0; j < n; j++)
		{
			mat1[i][j] = mat2[i][j] = (double)(rand() % n);
		}
	}

	cout << "Input:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	cout << "==================== Original method (literature):" << endl;
	cout << endl;

	Daub_StandardDecomposition(mat1, n, n, true);

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	Daub_StandardComposition(mat1, n, n, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	cout << "==================== Developed method:" << endl;
	cout << endl;

	cout << "Input:" << endl;
	cout << endl;
	printMatrix(mat2, n);
	cout << endl;

	Daub_StandardDecomposition(mat2, n, n, false, true);

	cout << "Optimal Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat2, n);
	cout << endl;

	//Daub_StandardNormalization(mat2, n);

	//cout << "Normalization Step Decomposition:" << endl;
	//cout << endl;
	//printMatrix(mat2, n);
	//cout << endl;

	Daub_StandardComposition(mat2, n, n, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat2, n);
	cout << endl;

	for (uint i = 0; i < n; i++)
	{
		delete[] mat1[i];
		delete[] mat2[i];
	}

	delete[] mat1;
	delete[] mat2;
}

void fundamentalTest9()
{
	cout << "********************" << endl;
	cout << "*     This fundamental (9) test is about the non standard normalized Daubechies" << endl;
	cout << "* (db2) decomposition and composition procedures. It uses a 16x16 matrix and" << endl;
	cout << "* performs no normalizing optimizations due to the fact it does not work on this" << endl;
	cout << "* approach." << endl;
	cout << "********************" << endl;
	cout << endl;

	uint n = 16;

	double **mat1 = new double*[n];

	for (uint i = 0; i < n; i++)
	{
		mat1[i] = new double[n];

		for (uint j = 0; j < n; j++)
		{
			mat1[i][j] = (double)(rand() % n);
		}
	}

	cout << "Input:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	cout << "==================== Developed method:" << endl;
	cout << endl;

	Daub_NonStandardDecomposition(mat1, n, n, true);

	cout << "Normalized Decomposition:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	Daub_NonStandardComposition(mat1, n, n, true);

	cout << "Normalized Composition:" << endl;
	cout << endl;
	printMatrix(mat1, n);
	cout << endl;

	for (uint i = 0; i < n; i++)
	{
		delete[] mat1[i];
	}

	delete[] mat1;
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
		case 2:
			fundamentalTest2();
			break;
		case 3:
			fundamentalTest3();
			break;
		case 4:
			fundamentalTest4();
			break;
		case 5:
			fundamentalTest5();
			break;
		case 6:
			fundamentalTest6();
			break;
		case 7:
			fundamentalTest7();
			break;
		case 8:
			fundamentalTest8();
			break;
		case 9:
			fundamentalTest9();
			break;
		default:
			cout << "\n\tFundamental test " << n << " not found!\n\n";
	}
}

NewArgs genNewArgs(int argc, char **argv, int takeOut)
{
	NewArgs result;

	result.argc = argc - takeOut;
	result.argv = (char**)malloc(sizeof(char**) * result.argc);

	for (int i = 0; i < result.argc; ++i)
	{
		result.argv[i] = argv[takeOut + i];
	}

	return result;
}

void freeNewArgs(NewArgs& newArgs)
{
	free(newArgs.argv);
}

void listAllTests()
{
	test12Desc();
	test13Desc();
	test14Desc();
	test15Desc();

	cout << endl;
}

int testIndexer(int argc, char **argv)
{
	int code = 0;

	if (argc < 2)
	{
		cout << endl;
		cout << "	* There is no default procedure for testing, you may " << endl;
		cout << "	like to list all tests by typing: ./tests.exe -l" << endl;
		cout << endl;
		return 1;
	}

	bool list = false;
	bool test = false;
	bool fundTest = false;

	if (strcmp(argv[1], "-l") == 0) list = true;
	else if (strcmp(argv[1], "-t") == 0) test = true;
	else if (strcmp(argv[1], "-ft") == 0) fundTest = true;

	if (list)
	{
		listAllTests();
	}
	else if (test)
	{
		if (argc < 3)
		{
			cout << endl;
			cout << "	* You have to tell the test number by typing: " << endl;
			cout << "	./tests.exe -t <test>" << endl;
			cout << endl;
			cout << "	* You can also list all possible tests by typing: " << endl;
			cout << "	./tests.exe -l" << endl;
			cout << endl;
			code = 1;
		}
		else
		{
			int ntest = atoi(argv[2]);
			NewArgs newArgs = genNewArgs(argc, argv, 3);

			switch(ntest)
			{
				case 1:
					test12(newArgs.argc, newArgs.argv);
					break;
				case 2:
					test13(newArgs.argc, newArgs.argv);
					break;
				case 3:
					test14(newArgs.argc, newArgs.argv);
					break;
				case 4:
					test15(newArgs.argc, newArgs.argv);
					break;
				default:
					cout << endl;
					cout << "	There is no test " << ntest << '.' << endl;
					cout << endl;
			}

			freeNewArgs(newArgs);
		}
	}
	else if (fundTest)
	{
		if (argc < 3)
		{
			cout << endl;
			cout << "	* You have to tell the test number by typing: " << endl;
			cout << "	./tests.exe -ft <test>" << endl;
			cout << endl;
			cout << "	* You can also list all possible tests by typing: " << endl;
			cout << "	./tests.exe -l" << endl;
			cout << endl;
			code = 1;
		}
		else
		{
			uint ntest = atoi(argv[2]);

			fundamentalTest(ntest);
		}
	}

	return code;
}
