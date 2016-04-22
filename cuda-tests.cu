#include "cuda-tests.h"

#ifdef __NVCC__

int test16()
{
	int n = 4;
	double *h_vec = new double[n];

	h_vec[0] = 9.0;
	h_vec[1] = 7.0;
	h_vec[2] = 5.0;
	h_vec[3] = 3.0;

	CUDA_Haar_Decomp(h_vec, n, false);

    cout << h_vec[0] << endl;
    cout << h_vec[1] << endl;
    cout << h_vec[2] << endl;
    cout << h_vec[3] << endl;

    delete [] h_vec;

	return 0;
}

#endif
