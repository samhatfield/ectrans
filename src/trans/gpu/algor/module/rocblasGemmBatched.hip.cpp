//
// Wrappers for rocblas_dgemm_batched and rocblas_sgemm_batched. 
//

#include <iostream>
#include <stdio.h>
#include "hip/hip_runtime_api.h"
#include "rocblas.h"

using namespace std;

#define ROCBLAS_CHECK(e)                                                \
  {                                                                     \
    rocblas_status err = (e);                                           \
    if (err != rocblas_status_success) {                                \
      fprintf(stderr, "ROCBLAS error: %s, line %d, %s: %i\n", __FILE__, \
              __LINE__, #e, err);                                       \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  }

bool alreadyAllocated_dgemm = false;
bool alreadyAllocated_sgemm = false;

double **d_Aarray_dgemm;
double **d_Barray_dgemm;
double **d_Carray_dgemm;

double **Aarray_dgemm;
double **Barray_dgemm;
double **Carray_dgemm;

float **d_Aarray_sgemm;
float **d_Barray_sgemm;
float **d_Carray_sgemm;

float **Aarray_sgemm;
float **Barray_sgemm;
float **Carray_sgemm;

rocblas_handle get_rocblas_handle() {
  static rocblas_handle handle;
  if (!handle) ROCBLAS_CHECK(rocblas_create_handle(&handle));
  return handle;
}

extern "C" void cublasDgemmBatched_wrapper (char transa, char transb,
                                            int m, int n, int k,
                                            double alpha,
                                            const double *A, int lda, int tda,
                                            const double *B, int ldb, int tdb,
                                            double beta, 
                                            double *C, int ldc, int tdc,
                                            int batchCount)
{
  // Determine whether either input array should be transposed
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  if (transa == 'T' || transa == 't') op_t1 = rocblas_operation_transpose;
  if (transb == 'T' || transb == 't') op_t2 = rocblas_operation_transpose;

  rocblas_handle handle = get_rocblas_handle();

  if (!alreadyAllocated_dgemm){
    hipHostMalloc(&Aarray_dgemm, batchCount*sizeof(double*));
    hipHostMalloc(&Barray_dgemm, batchCount*sizeof(double*));
    hipHostMalloc(&Carray_dgemm, batchCount*sizeof(double*));

    hipMalloc(&d_Aarray_dgemm, batchCount*sizeof(double*));
    hipMalloc(&d_Barray_dgemm, batchCount*sizeof(double*));
    hipMalloc(&d_Carray_dgemm, batchCount*sizeof(double*));
   }
  alreadyAllocated_dgemm = true;

  for (int i=0; i < batchCount; i++) {
    Aarray_dgemm[i] = (double*) &(A[i*lda*tda]);
    Barray_dgemm[i] = (double*) &(B[i*ldb*tdb]);
    Carray_dgemm[i] = (double*) &(C[i*ldc*tdc]);
  }

  hipMemcpy(d_Aarray_dgemm, Aarray_dgemm, batchCount*sizeof(double*), hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_dgemm, Barray_dgemm, batchCount*sizeof(double*), hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_dgemm, Carray_dgemm, batchCount*sizeof(double*), hipMemcpyHostToDevice);

  ROCBLAS_CHECK(rocblas_dgemm_batched(handle, op_t1, op_t2, m, n, k,
                                      &alpha,
                                      (const double**) d_Aarray_dgemm, lda,
                                      (const double**) d_Barray_dgemm, ldb,
                                      &beta,
                                      (double**) d_Carray_dgemm, ldc,
                                      batchCount));
  hipDeviceSynchronize();
}

extern "C" void cublasSgemmBatched_wrapper(char transa, char transb,
                                           int m, int n, int k,
                                           float alpha,
                                           const float *A, int lda, int tda,
                                           const float *B, int ldb, int tdb,
                                           float beta,
                                           float *C, int ldc, int tdc,
                                           int batchCount)
{
  // Determine whether either input array should be transposed
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  if (transa == 'T' || transa == 't') op_t1 = rocblas_operation_transpose;
  if (transb == 'T' || transb == 't') op_t2 = rocblas_operation_transpose;

  rocblas_handle handle = get_rocblas_handle();

  if (!alreadyAllocated_sgemm) {
    hipHostMalloc(&Aarray_sgemm, batchCount*sizeof(float*));
    hipHostMalloc(&Barray_sgemm, batchCount*sizeof(float*));
    hipHostMalloc(&Carray_sgemm, batchCount*sizeof(float*));

    hipMalloc(&d_Aarray_sgemm, batchCount*sizeof(float*));
    hipMalloc(&d_Barray_sgemm, batchCount*sizeof(float*));
    hipMalloc(&d_Carray_sgemm, batchCount*sizeof(float*));
  }
  alreadyAllocated_sgemm = true;

  for (int i = 0; i < batchCount; i++) {
    Aarray_sgemm[i] = (float*) &(A[i*lda*tda]);
    Barray_sgemm[i] = (float*) &(B[i*ldb*tdb]);
    Carray_sgemm[i] = (float*) &(C[i*ldc*tdc]);
  }

  hipMemcpy(d_Aarray_sgemm, Aarray_sgemm, batchCount*sizeof(float*), hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_sgemm, Barray_sgemm, batchCount*sizeof(float*), hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_sgemm, Carray_sgemm, batchCount*sizeof(float*), hipMemcpyHostToDevice);

  ROCBLAS_CHECK(rocblas_sgemm_batched(handle, op_t1, op_t2, m, n, k,
                                      &alpha,
                                      (const float**) d_Aarray_sgemm, lda,
                                      (const float**) d_Barray_sgemm, ldb,
                                      &beta,
                                      (float**) d_Carray_sgemm, ldc,
                                      batchCount));

  hipDeviceSynchronize();
}

// This should be implemented eventually, and also for DGEMM
// extern "C" void cublasSgemmStridedBatched_wrapper(char transa, char transb,
//                                                   int m, int n, int k,
//                                                   float alpha,
//                                                   const float *A, int lda, int tda,
//                                                   const float *B, int ldb, int tdb,
//                                                   float beta,
//                                                   float *C, int ldc, int tdc,
//                                                   int batchCount)
// {
//   // Determine whether either input array should be transposed
//   rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
//   if (transa=='T' || transa=='t') op_t1 = rocblas_operation_transpose;
//   if (transb=='T' || transb=='t') op_t2 = rocblas_operation_transpose;

//   rocblas_handle handle = get_rocblas_handle();

//   ROCBLAS_CHECK(rocblas_sgemm_strided_batched(handle, op_t1, op_t2, m, n, k,
//                                               &alpha,
//                                               (const float *) A, lda, tda,
//                                               (const float *) B, ldb, tdb,
//                                               &beta,
//                                               (float*) C, ldc, tdc,
//                                               batchCount));
// }