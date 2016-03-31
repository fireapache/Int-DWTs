#include "cuda-tests.h"

#ifdef __NVCC__

int test16()
{
	cudaError_t err = cudaSuccess;
	int n = 1024;
	int size = sizeof(double) * n;

	double *h_vec = new double[n];
	double *d_vec = NULL;

	err = cudaMalloc((void**)&d_vec, size);

	if (err != cudaSuccess)
    {
        cout << "Failed to allocate device vector!" << endl;
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(d_vec, h_vec, size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        cout << "Failed to copy vector A from host to device!" << endl;
        exit(EXIT_FAILURE);
    }

	int threads = 256;
	int blocks = n / threads;

	CUDA_Haar_Decomp <<< threads, blocks >>> (d_vec, n, true);

	cudaDeviceSynchronize();

	err = cudaMemcpy(h_vec, d_vec, size, cudaMemcpyDeviceToHost);

	if (err != cudaSuccess)
    {
        cout << "Failed to copy vector C from device to host!" << endl;
        exit(EXIT_FAILURE);
    }

    cudaFree(d_vec);

    cout << h_vec[0] << endl;

    delete [] h_vec;

    err = cudaDeviceReset();

    if (err != cudaSuccess)
    {
        cout << "Failed to deinitialize the device!" << endl;
        exit(EXIT_FAILURE);
    }

	return 0;
}

#endif
