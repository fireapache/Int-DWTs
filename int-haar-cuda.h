#pragma once

#if __NVCC__

#include <cuda_runtime.h>
//#include <interval.hpp>
#include <iostream>

//using namespace cxsc;
using namespace std;

void CUDA_Haar_Decomp(double *h_vec, int n, bool normal);

#endif