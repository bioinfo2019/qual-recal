/*
 * nrmap.h
 *
 *  Created on: Mar 11, 2017
 *      Author: eliot
 */

#ifndef NRMAP_H_
#define NRMAP_H_


#endif /* NRMAP_H_ */

#include <string>
#include <zlib.h>
#include <pthread.h>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <cstring>
#include <algorithm>

#include <numeric>
#include "bamlib.hpp"


#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <assert.h>

#include <cstdio>
#include <cstdlib>


float* getInvariants(unsigned char* h_readSequences, int rows, int cols);
void gpuMatMult(int n, int k, float* B, int ldb, float* h_lsq, int ldc, int nReads, int readLength);
int getCorrelationCoefficients(int rows, int cols, float* h_corrCoefs, unsigned char* h_readQuals);
