#pragma once

#if __NVCC__

#include <cuda_runtime.h>

void __global__ CUDA_Haar_Decomp(double *vec, int n, bool normal);
void __device__ CUDA_Haar_DecompStep(double *vec, int n, bool normal);

#endif