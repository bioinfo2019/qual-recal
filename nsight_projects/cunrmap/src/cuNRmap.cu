/*
 ============================================================================
 Name        : cuNRmap.cu
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include "nrmap.h"
//#include "jacobi_eigenvalue.hpp"

static void CheckCudaErrorAux (const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

#define THREADBLOCK_SIZE 64
#define ARRAY_SIZE 768
#define ROWS 20
#define COLS 10

#define nTPB 64

#define ROW_MAJOR 0
#define COL_MAJOR 1
#define PSEUDO_INVERSE_MATRIX_LEN 200;
#define PSEUDO_INVERSE_MATRIX_LDA 2;


void gpuMatMult(int n, int k, float* B, int ldb, float* h_lsq, int ldc, int nReads, int readLength)
{
	int m = 2;

	float pseudo_inv_100bp[200] = {0.04, -0.0005940594, 0.0393939394, -0.0005820582, 0.0387878788, -0.000570057, 0.0381818182, -0.0005580558, 0.0375757576, -0.0005460546, 0.036969697, -0.0005340534, 0.0363636364, -0.0005220522, 0.0357575758, -0.000510051, 0.0351515152, -0.0004980498, 0.0345454545, -0.0004860486, 0.0339393939, -0.0004740474, 0.0333333333, -0.0004620462, 0.0327272727, -0.000450045, 0.0321212121, -0.0004380438, 0.0315151515, -0.0004260426, 0.0309090909, -0.0004140414, 0.0303030303, -0.0004020402, 0.0296969697, -0.000390039, 0.0290909091, -0.0003780378, 0.0284848485, -0.0003660366, 0.0278787879, -0.0003540354, 0.0272727273, -0.0003420342, 0.0266666667, -0.000330033, 0.0260606061, -0.0003180318, 0.0254545455, -0.0003060306, 0.0248484848, -0.0002940294, 0.0242424242, -0.0002820282, 0.0236363636, -0.000270027, 0.023030303, -0.0002580258, 0.0224242424, -0.0002460246, 0.0218181818, -0.0002340234, 0.0212121212, -0.0002220222, 0.0206060606, -0.000210021, 0.02, -0.0001980198, 0.0193939394, -0.0001860186, 0.0187878788, -0.0001740174, 0.0181818182, -0.0001620162, 0.0175757576, -0.000150015, 0.016969697, -0.0001380138, 0.0163636364, -0.0001260126, 0.0157575758, -0.0001140114, 0.0151515152, -0.0001020102, 0.0145454545, -0.000090009, 0.0139393939, -7.80078007800747E-05, 0.0133333333, -6.6006600660063E-05, 0.0127272727, -5.40054005400512E-05, 0.0121212121, -4.20042004200394E-05, 0.0115151515, -0.000030003, 0.0109090909, -1.80018001800158E-05, 0.0103030303, -6.000600060004E-06, 0.0096969697, 6.0006000600078E-06, 0.0090909091, 1.80018001800196E-05, 0.0084848485, 0.000030003, 0.0078787879, 4.20042004200432E-05, 0.0072727273, 5.4005400540055E-05, 0.0066666667, 6.60066006600668E-05, 0.0060606061, 7.80078007800786E-05, 0.0054545455, 0.000090009, 0.0048484848, 0.0001020102, 0.0042424242, 0.0001140114, 0.0036363636, 0.0001260126, 0.003030303, 0.0001380138, 0.0024242424, 0.000150015, 0.0018181818, 0.0001620162, 0.0012121212, 0.0001740174, 0.0006060606, 0.0001860186, 1.17961196366423E-16, 0.0001980198, -0.0006060606, 0.000210021, -0.0012121212, 0.0002220222, -0.0018181818, 0.0002340234, -0.0024242424, 0.0002460246, -0.003030303, 0.0002580258, -0.0036363636, 0.000270027, -0.0042424242, 0.0002820282, -0.0048484848, 0.0002940294, -0.0054545455, 0.0003060306, -0.0060606061, 0.0003180318, -0.0066666667, 0.000330033, -0.0072727273, 0.0003420342, -0.0078787879, 0.0003540354, -0.0084848485, 0.0003660366, -0.0090909091, 0.0003780378, -0.0096969697, 0.000390039, -0.0103030303, 0.0004020402, -0.0109090909, 0.0004140414, -0.0115151515, 0.0004260426, -0.0121212121, 0.0004380438, -0.0127272727, 0.000450045, -0.0133333333, 0.0004620462, -0.0139393939, 0.0004740474, -0.0145454545, 0.0004860486, -0.0151515152, 0.0004980498, -0.0157575758, 0.000510051, -0.0163636364, 0.0005220522, -0.016969697, 0.0005340534, -0.0175757576, 0.0005460546, -0.0181818182, 0.0005580558, -0.0187878788, 0.000570057, -0.0193939394, 0.0005820582, -0.02, 0.0005940594};

	float *d_psinv, *d_readquals, *d_lsq;

	float alpha = 1.0;
	float beta = 0.0;

	cudaMalloc(&d_psinv, sizeof(float)*200);
	cudaMemcpy(d_psinv, pseudo_inv_100bp, sizeof(float)*200, cudaMemcpyHostToDevice);

	cudaMalloc(&d_readquals, sizeof(float)*nReads*readLength);
	cudaMemcpy(d_readquals, B, sizeof(float)*nReads*readLength, cudaMemcpyHostToDevice);

	cudaMalloc(&d_lsq, sizeof(float)*2*nReads);
	cudaMemset(d_lsq, 0, sizeof(float)*2*nReads);

	cublasHandle_t myhandle;
	cublasStatus_t cublas_result;

    cublas_result = cublasCreate(&myhandle);
    assert(cublas_result == CUBLAS_STATUS_SUCCESS);

    cublas_result = cublasSgemmStridedBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, d_psinv, 2, 0, d_readquals, readLength, readLength, &beta, d_lsq, ldc, 2, nReads);

	cudaMemcpy(h_lsq, d_lsq, sizeof(float)*2*nReads, cudaMemcpyDeviceToHost);

	cudaFree(d_psinv);

	cudaFree(d_readquals);

	cudaFree(d_lsq);

	assert(cublas_result == CUBLAS_STATUS_SUCCESS);

}


template <int select, typename T>
__global__ void vec_mat_row_add(const unsigned int height, const unsigned int width, T* matrix, const T* vector)
{
    // get the current element index for the thread
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < height*width)
    {
        // sum the current element with the
    if (select == ROW_MAJOR)
        matrix[idx] += vector[idx%width];
    else // COL_MAJOR
        matrix[idx] += vector[idx/height];
    }
}

/**
 * CUDA kernel that computes reciprocal values for a given vector
 */
__global__ void reciprocalKernel(float *data, unsigned vectorSize) {
	unsigned idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx < vectorSize)
		data[idx] = 1.0/data[idx];
}

/**
 * Host function that copies the data and launches the work on GPU
 */
float *gpuReciprocal(float *data, unsigned size)
{
	float *rc = new float[size];
	float *gpuData;

	CUDA_CHECK_RETURN(cudaMalloc((void **)&gpuData, sizeof(float)*size));
	CUDA_CHECK_RETURN(cudaMemcpy(gpuData, data, sizeof(float)*size, cudaMemcpyHostToDevice));
	
	static const int BLOCK_SIZE = 256;
	const int blockCount = (size+BLOCK_SIZE-1)/BLOCK_SIZE;
	reciprocalKernel<<<blockCount, BLOCK_SIZE>>> (gpuData, size);

	CUDA_CHECK_RETURN(cudaMemcpy(rc, gpuData, sizeof(float)*size, cudaMemcpyDeviceToHost));
	CUDA_CHECK_RETURN(cudaFree(gpuData));
	return rc;
}

float *cpuReciprocal(float *data, unsigned size)
{
	float *rc = new float[size];
	for (unsigned cnt = 0; cnt < size; ++cnt) rc[cnt] = 1.0/data[cnt];
	return rc;
}


void initialize(float *data, unsigned size)
{
	for (unsigned i = 0; i < size; ++i)
		data[i] = .5*(i+1);
}


__device__ __forceinline__ int no_bank_conflict_index(int thread_id,
                                                      int logical_index)
{
    return logical_index * THREADBLOCK_SIZE + thread_id;
}


__device__  __forceinline__  void r8mat_diag_get_vector ( int n, float a[], float v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N*N], the N by N matrix.
//
//    Output, double V[N], the diagonal entries
//    of the matrix.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return;
}
//****************************************************************************80

__device__  __forceinline__  void r8mat_identity ( int n, float a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IDENTITY sets the square matrix A to the identity.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double A[N*N], the N by N identity matrix.
//
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
//****************************************************************************80

__device__  __forceinline__  void jacobi_eigenvalue (float a[], int it_max, float v[],
  float d[], int &it_num, int &rot_num, int alphabet_size)

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
//
//  Discussion:
//
//    This function computes the eigenvalues and eigenvectors of a
//    real symmetric matrix, using Rutishauser's modfications of the classical
//    Jacobi rotation method with threshold pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2013
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix, which must be square, real,
//    and symmetric.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, double V[N*N], the matrix of eigenvectors.
//
//    Output, double D[N], the eigenvalues, in descending order.
//
//    Output, int &IT_NUM, the total number of iterations.
//
//    Output, int &ROT_NUM, the total number of rotations.
//
{
	float bw[4];
	float c;
	float g;
	float gapq;
	float h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  float s;
  float t;
  float tau;
  float term;
  float termp;
  float termq;
  float theta;
  float thresh;
  float w;
  float zw[4];

  int n = 4;

  r8mat_identity ( n, v );

  r8mat_diag_get_vector ( n, a, d );

  for ( i = 0; i < n; i++ )
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;

  while ( it_num < it_max )
  {
    it_num = it_num + 1;
//
//  The convergence threshold is based on the size of the elements in
//  the strict upper triangle of the matrix.
//
    thresh = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < j; i++ )
      {
        thresh = thresh + a[i+j*n] * a[i+j*n];
      }
    }

    thresh = sqrt ( thresh ) / ( double ) ( 4 * n );

    if ( thresh == 0.0 )
    {
      break;
    }

    for ( p = 0; p < n; p++ )
    {
      for ( q = p + 1; q < n; q++ )
      {
        gapq = 10.0 * fabs ( a[p+q*n] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
//
//  Annihilate tiny offdiagonal elements.
//
        if ( 4 < it_num &&
             termp == fabs ( d[p] ) &&
             termq == fabs ( d[q] ) )
        {
          a[p+q*n] = 0.0;
        }
//
//  Otherwise, apply a rotation.
//
        else if ( thresh <= fabs ( a[p+q*n] ) )
        {
          h = d[q] - d[p];
          term = fabs ( h ) + gapq;

          if ( term == fabs ( h ) )
          {
            t = a[p+q*n] / h;
          }
          else
          {
            theta = 0.5 * h / a[p+q*n];
            t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 )
            {
              t = - t;
            }
          }
          c = 1.0 / sqrt ( 1.0 + t * t );
          s = t * c;
          tau = s / ( 1.0 + c );
          h = t * a[p+q*n];
//
//  Accumulate corrections to diagonal elements.
//
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[p+q*n] = 0.0;
//
//  Rotate, using information from the upper triangle of A only.
//
          for ( j = 0; j < p; j++ )
          {
            g = a[j+p*n];
            h = a[j+q*n];
            a[j+p*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = p + 1; j < q; j++ )
          {
            g = a[p+j*n];
            h = a[j+q*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = q + 1; j < n; j++ )
          {
            g = a[p+j*n];
            h = a[q+j*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[q+j*n] = h + s * ( g - h * tau );
          }
//
//  Accumulate information in the eigenvector matrix.
//
          for ( j = 0; j < n; j++ )
          {
            g = v[j+p*n];
            h = v[j+q*n];
            v[j+p*n] = g - s * ( h + g * tau );
            v[j+q*n] = h + s * ( g - h * tau );
          }
          rot_num = rot_num + 1;
        }
      }
    }

    for ( i = 0; i < n; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
//
//  Restore upper triangle of input matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }
//
//  Ascending sort the eigenvalues and eigenvectors.
//
  for ( k = 0; k < alphabet_size - 1; k++ )
  {
    m = k;
    for ( l = k + 1; l < alphabet_size; l++ )
    {
      if ( d[l] < d[m] )
      {
        m = l;
      }
    }

    if ( m != k )
    {
      t    = d[m];
      d[m] = d[k];
      d[k] = t;
      //for ( i = 0; i < n; i++ )
      //{
      //  w        = v[i+m*n];
      //  v[i+m*n] = v[i+k*n];
      //  v[i+k*n] = w;
      //}
    }
  }

  return;
}
//****************************************************************************80

   // assumes matrix indices start from 0 (0,1,2 and 3)
__device__ float determinant(float m[]) {


	return
			 m[3] * m[6] * m[9] * m[12] - m[2] * m[7] * m[9] * m[12] -
			 m[3] * m[5] * m[10] * m[12] + m[1] * m[7] * m[10] * m[12] +
			 m[2] * m[5] * m[11] * m[12] - m[1] * m[6] * m[11] * m[12] -
			 m[3] * m[6] * m[8] * m[13] + m[2] * m[7] * m[8] * m[13] +
			 m[3] * m[4] * m[10] * m[13] - m[0] * m[7] * m[10] * m[13] -
			 m[2] * m[4] * m[11] * m[13] + m[0] * m[6] * m[11] * m[13] +
			 m[3] * m[5] * m[8] * m[14] - m[1] * m[7] * m[8] * m[14] -
			 m[3] * m[4] * m[9] * m[14] + m[0] * m[7] * m[9] * m[14] +
			 m[1] * m[4] * m[11] * m[14] - m[0] * m[5] * m[11] * m[14] -
			 m[2] * m[5] * m[8] * m[15] + m[1] * m[6] * m[8] * m[15] +
			 m[2] * m[4] * m[9] * m[15] - m[0] * m[6] * m[9] * m[15] -
			 m[1] * m[4] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];

}

__device__ float trace(float m[]) {

	return m[0] + m[5] + m[10] + m[15];

}

__global__ void corrcoef_gpu(int height, int width, float* matrix, const unsigned char* qual_scores)
{
	__shared__ unsigned char shm_qual_scores[THREADBLOCK_SIZE*ARRAY_SIZE];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid < height)
	{
		unsigned char *read_qual_scores = (unsigned char*)((unsigned char*)qual_scores + (tid)*width);

		for (unsigned int i = 0;i < width;++i)
		{
			shm_qual_scores[no_bank_conflict_index(threadIdx.x, i)] = read_qual_scores[i];
		}

		float sumOfY = 0;
		float sumOfXY = 0;
		float sumOfX = 0;
		float meanOfX = 0;
		float meanOfY = 0;
		float sumOfXSquared = 0;
		float sumOfYSquared = 0;

		for (unsigned int i = 0;i < width;++i)
		{
			sumOfY += shm_qual_scores[no_bank_conflict_index(threadIdx.x, i)];
			sumOfXY += (i + 1) * shm_qual_scores[no_bank_conflict_index(threadIdx.x, i)];
			sumOfX += (i + 1);
			sumOfXSquared += (i+1)*(i+1);
			sumOfYSquared += shm_qual_scores[no_bank_conflict_index(threadIdx.x, i)]*shm_qual_scores[no_bank_conflict_index(threadIdx.x, i)];

		}

		meanOfY = (float)sumOfY/(float)width;
		meanOfX = (float)sumOfX / (float)width;


		matrix[tid] = (sumOfXY - width*meanOfX*meanOfY)/(sqrt(sumOfXSquared - width*meanOfX*meanOfX)*sqrt(sumOfYSquared - width*meanOfY*meanOfY));

	}

}

int getCorrelationCoefficients(int rows, int cols, float* h_corrCoefs, unsigned char* h_readQuals)
{

	//Fix the block size to an appropriate value
	dim3 block(THREADBLOCK_SIZE);
	dim3 grid;
	grid.x = (rows + block.x - 1)/block.x;

	unsigned char *d_readQuals;
	float *d_corrCoefs;

	cudaMalloc(&d_corrCoefs, rows*sizeof(float));
	cudaMemset(d_corrCoefs, 0, rows*sizeof(float));

	cudaMalloc(&d_readQuals, sizeof(unsigned char)*cols*rows);

	cudaMemcpy(d_readQuals, h_readQuals, sizeof(unsigned char)*cols*rows, cudaMemcpyHostToDevice);


	corrcoef_gpu<<<grid.x, block.x>>>(rows, cols, d_corrCoefs, d_readQuals);

	cudaDeviceSynchronize();
	cudaMemcpy(h_corrCoefs, d_corrCoefs, rows*sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(d_readQuals);
	cudaFree(d_corrCoefs);

	return 0;

}

__global__ void calc_invariants_gpu(int height, int width, float* matrix, const unsigned char* dna_seqs)
{

	__shared__ unsigned char idxs[THREADBLOCK_SIZE*ARRAY_SIZE];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	const int max_iter = 1000;

	//unsigned char idxs[256] = {0};

	if (tid < height)
	{
		//start out with the condensed matrix initialized to the 4x4 identity matrix
		float cm[16*sizeof(float)] = {0};

		cm[0] = 1;
		cm[5] = 1;
		cm[10] = 1;
		cm[15] = 1;

		unsigned char alphabet[4] = {0};
		alphabet[0] = 0;
		alphabet[1] = 0;
		alphabet[2] = 0;
		alphabet[3] = 0;

		// update the pointer to point to the beginning of the row
	    unsigned char *dna_seq = (unsigned char*)((unsigned char*)dna_seqs + (tid)*width);

		alphabet[0] = dna_seq[0];

		int alphabet_size = 0;

		float gcCount = 0;
		float sequenceEffectiveLength = 0;

		for (unsigned int i = 0;i < width;++i)
		{
			if(dna_seq[i] != 'N')
				sequenceEffectiveLength += 1;

			if(dna_seq[i] == 'C' || dna_seq[i] == 'G')
				gcCount += 1.0;
		}


		float gcContent = gcCount/(sequenceEffectiveLength);

		for(int l = 1; l < width; ++l)
		{
			if(dna_seq[l] == (unsigned char) 'N')
				continue;

			if(dna_seq[l] != alphabet[0] && dna_seq[l] != alphabet[1] && dna_seq[l] != alphabet[2] && dna_seq[l] != alphabet[3]) //need to cover the case where a residue is N
			{
				alphabet[alphabet_size] = dna_seq[l];
				alphabet_size += 1;
			}
		}


		int base1_counter = 0;

		//loop through letters - outer pass
		int base1 = 0;
		for (base1 = 0;base1 < alphabet_size; ++base1)
		{
			base1_counter = 0;
			//collect indexes of first letter and store in shared memory
			for (int i = 0;i < width; ++i)
			{
				if(dna_seq[i] == (unsigned char) 'N')
					continue;

				if (alphabet[base1] == dna_seq[i])
				{
					idxs[no_bank_conflict_index(threadIdx.x, base1_counter)] = i;
					base1_counter += 1;
				}
			}


			int base2 = 0;
			//loop through letters - inner pass
			for (base2 = 0;base2 < alphabet_size; ++base2)
			{
				int base2_counter = base1_counter;
				//collect indexes of second letter and store in shared memory
				for (int j = 0;j < width; ++j)
				{
					if(dna_seq[j] == (unsigned char) 'N')
						continue;

					if (alphabet[base2] == dna_seq[j])
					{
						idxs[no_bank_conflict_index(threadIdx.x, base2_counter)] = j;

						base2_counter += 1;
					}
				}

				//this block sums up the distances in the distance matrix for each letter pair
				//note that since we only want the average distance from one type of nucleotide
				//to another, we don't need to actually create the matrix. Just accumulate the
				//sum of the distances and take the mean. This saves shit tons of memory!
				float dist = 0;

				for(int row = 0;row < base1_counter; ++row)
					for(int col = base1_counter;col < base2_counter;++col)
						dist += fabs((float)idxs[no_bank_conflict_index(threadIdx.x, row)] - (float)idxs[no_bank_conflict_index(threadIdx.x, col)]);

				cm[base2*4 + base1] = dist/(((float)base1_counter)*((float)(base2_counter-base1_counter)));
				//matrix[base2*alphabet_size + base1] = dist/(((float)base1_counter)*((float)(base2_counter-base1_counter)));

			}

		}


		float mx = 0;
		float mn = 100000;
		int idx = 0;

		for(int row = 0;row < alphabet_size;++row)
			for(int col = 0;col < alphabet_size;++col)
		{
			idx = row * 4 + col;

			if (cm[idx] > mx)
				mx = cm[idx];
			if (cm[idx] < mn)
				mn = cm[idx];
		}

		int count = 0;

		for(int row = 0;row < alphabet_size;++row)
			for(int col = 0;col < alphabet_size;++col)
		{
			idx = row * 4 + col;
			cm[idx] = (cm[idx] - mn)/(mx - mn);
		}

		int it_num = 0;
		int rot_num = 0;
		float v[16*sizeof(float)] = {0};
		float d[4*sizeof(float)] = {0};
		jacobi_eigenvalue(cm, max_iter, v, d, it_num, rot_num, alphabet_size);
		matrix[tid * INVARIANTS_BUFFER_WIDTH] = d[alphabet_size - 1];
		matrix[tid * INVARIANTS_BUFFER_WIDTH + 1] = d[alphabet_size - 2];

		// handle the case where our condensed matrix has only the upper left 2x2 or 3x3 elements filled in
		// otherwise the 1's in the lower part of the diagonal will inflate the value of the trace.
		// Those 1's prevent the value of the determinant for the 2x2 and 3x3 cases from being zero

		matrix[tid * INVARIANTS_BUFFER_WIDTH + 2] = trace(cm) - (4 - alphabet_size);

		matrix[tid * INVARIANTS_BUFFER_WIDTH + 3] = d[0] * d[1] * d[2] * d[3]; //determinant - product of eigenvalues
		matrix[tid * INVARIANTS_BUFFER_WIDTH + 4] = gcContent;


	}

	//__syncthreads();

}


float* getInvariants(unsigned char* h_readSequences, int rows, int cols)
{

	//Fix the block size to an appropriate value
	dim3 block(THREADBLOCK_SIZE);
	dim3 grid;
	grid.x = (rows + block.x - 1)/block.x;

	unsigned char *d_readSequences;
	float *d_invariants, *h_invariants;

	const int invsz = rows*INVARIANTS_BUFFER_WIDTH*sizeof(float);

	h_invariants = (float *)malloc(invsz);
	cudaMalloc(&d_invariants, invsz);
	cudaMemset(d_invariants, 0, invsz);

	cudaMalloc(&d_readSequences, sizeof(unsigned char)*cols*rows);

	cudaMemcpy(d_readSequences, h_readSequences, sizeof(unsigned char)*cols*rows, cudaMemcpyHostToDevice);

	calc_invariants_gpu<<<grid.x, block.x>>>(rows, cols, d_invariants, d_readSequences);

	cudaDeviceSynchronize();
	cudaMemcpy(h_invariants, d_invariants, invsz, cudaMemcpyDeviceToHost);

	cudaFree(d_readSequences);
	cudaFree(d_invariants);

	return h_invariants;

}


void calc_invariants_cpu(const unsigned int height, const unsigned int width, float* matrix, const unsigned* dna_seq)
{
	unsigned alphabet[4] = {0};
	//extern __shared__ unsigned idxs[];

	//for (int feat = 0;feat < height;++feat)

	//{
		unsigned idxs[92] = {0};

		//unsigned *dna_seq;

		//int a = feat*2; //???????

		///dna_seq = &dna_seqs[a];

		alphabet[0] = dna_seq[0];

		int alphabet_size = 0;

		for(int l = 1; l < width; ++l)
		{
			if(dna_seq[l] != alphabet[0] && dna_seq[l] != alphabet[1] && dna_seq[l] != alphabet[2] && dna_seq[l] != alphabet[3]) //need to cover the case where a residue is N
			{
				alphabet[alphabet_size] = dna_seq[l];
				alphabet_size += 1;
			}
		}

		unsigned base1_counter = 0;

		//loop through letters - outer pass
		unsigned base1 = 0;
		for (base1 = 0;base1 < alphabet_size; ++base1)
		{
			base1_counter = 0;
			//collect indexes of first letter and store in shared memory
			for (unsigned i = 0;i < width; ++i)
			{
				if (alphabet[base1] == dna_seq[i])
				{
					idxs[base1_counter] = i;
					base1_counter += 1;
				}
			}

			unsigned base2_counter = 0;
			unsigned base2 = 0;
			//loop through letters - inner pass
			for (base2 = 0;base2 < alphabet_size; ++base2)
			{
				base2_counter = 0;
				//collect indexes of second letter and store in shared memory
				for (unsigned j = 0;j < width; ++j)
				{
					if (alphabet[base2] == dna_seq[j])
					{
						idxs[base1_counter + base2_counter] = j;
						base2_counter += 1;
					}
				}

				//this block sums up the distance in the distance matrix for each letter pair
				//note that since we only want the average distance from one type of nucleotide
				//to another, we don't need to actually create the matrix. Just accumulate the
				//sum of the distances and take the mean. This saves shit tons of memory!
				float dist = 0;

				for(int row = 0;row < base1_counter;++row)
					for(int col = base1_counter;col < base1_counter + base2_counter;++col)
						dist += fabs((float)idxs[row] - (float)idxs[col]);

				matrix[base2*alphabet_size + base1] = dist/(((float)base1_counter)*((float)base2_counter));
			}

		}
	//}
}


/**
 * Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 */
static void CheckCudaErrorAux (const char *file, unsigned line, const char *statement, cudaError_t err)
{
	if (err == cudaSuccess)
		return;
	std::cerr << statement<<" returned " << cudaGetErrorString(err) << "("<<err<< ") at "<<file<<":"<<line << std::endl;
	exit (1);
}

