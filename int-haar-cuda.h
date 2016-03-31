#pragma once

#if __NVCC__

#include <cuda_runtime.h>

template <typename T>
__global__ void CUDA_Haar_Decomp(T *vec, int n, bool normal)
{
	vec[0] = T(69.5);
}

template <typename T>
__device__ void CUDA_Haar_DecompStep(T *vec, int n, bool normal)
{
	
}

#endif