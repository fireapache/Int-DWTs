#include "int-haar-cuda.h"

#ifdef __NVCC__

template <typename T>
__global__ void _CUDA_Haar_Normalization_Step(T *data, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < n)
	{
		data[i] = data[i] / sqrt(float(n));
	}
}

template <typename T>
__global__ void _CUDA_Haar_DecompStep(T *data, T *result, int n, T div)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < n)
	{
		result[i]     = (data[2*i] + data[2*i + 1]) / div;
		result[i + n] = (data[2*i] - data[2*i + 1]) / div;
	}
}

template <typename T>
__global__ void _CUDA_Haar_Next_Data_Level(T *data, T *result, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < n)
	{
		data[i] = result[i];
	}
}

template <typename T>
void _CUDA_Haar_Decomp(T *h_vec, int n, bool normal)
{
	cudaError_t err = cudaSuccess;
	T *d_vec = NULL;
	T *r_vec = NULL;
	int size = sizeof(T) * n;
	int threads = 256;
	int blocks;

	if (threads > n) threads = n;
	blocks = n / threads;

	err = cudaMalloc((void**)&d_vec, size);

	if (err != cudaSuccess)
	{
		cout << "Failed to allocate device vector" << endl;
		exit(EXIT_FAILURE);
	}

    err = cudaMemcpy(d_vec, h_vec, size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        cout << "Failed to copy vector A from host to device!" << endl;
        exit(EXIT_FAILURE);
    }

	if (normal)
	{
		_CUDA_Haar_Normalization_Step <<< blocks, threads >>> (d_vec, n);
		cudaDeviceSynchronize();
	}

	int n_step = n;
	int size_step;
	T div;

	if (normal) div = sqrt(T(2.0));
	else        div = T(2.0);

	T *temp = new T[n];

	while (n_step > 1)
	{
		size_step = sizeof(T) * n_step;
		err = cudaMalloc((void**)&r_vec, size_step);

		if (err != cudaSuccess)
		{
			cout << "Failed to allocate result device vector!" << endl;
			exit(EXIT_FAILURE);
		}

		if (threads > n_step / 2) threads = n_step / 2;
		blocks = (n_step / 2) / threads;

		cudaMemcpy(temp, d_vec, size, cudaMemcpyDeviceToHost);

		cout << temp[0] << endl;
	    cout << temp[1] << endl;
	    cout << temp[2] << endl;
	    cout << temp[3] << endl;
	    cout << "-----" << endl;

		_CUDA_Haar_DecompStep <<< blocks, threads >>> (d_vec, r_vec, n_step/2, div);
		cudaDeviceSynchronize();

		threads *= 2;

		if (threads > n) threads = n;

		_CUDA_Haar_Next_Data_Level <<< blocks, threads >>> (d_vec, r_vec, n_step);
		cudaDeviceSynchronize();

		cudaMemcpy(temp, d_vec, size, cudaMemcpyDeviceToHost);

		cout << temp[0] << endl;
	    cout << temp[1] << endl;
	    cout << temp[2] << endl;
	    cout << temp[3] << endl;
	    cout << "-----" << endl;

		err = cudaFree(r_vec);

		if (err != cudaSuccess)
		{
			cout << "Failed to free result vector on device!" << endl;
			exit(EXIT_FAILURE);
		}

		n_step /= 2;
	}

	delete [] temp;

	err = cudaMemcpy(h_vec, d_vec, size, cudaMemcpyDeviceToHost);

	if (err != cudaSuccess)
    {
        cout << "Failed to copy vector C from device to host!" << endl;
        exit(EXIT_FAILURE);
    }

    cudaFree(d_vec);

    err = cudaDeviceReset();

    if (err != cudaSuccess)
    {
        cout << "Failed to deinitialize the device!" << endl;
        exit(EXIT_FAILURE);
    }

}

void CUDA_Haar_Decomp(double *h_vec, int n, bool normal)
{
	_CUDA_Haar_Decomp(h_vec, n, normal);
}


#endif
