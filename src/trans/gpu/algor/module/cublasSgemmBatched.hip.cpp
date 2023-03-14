//
// Wrapper for hipblasSgemm function. 
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

bool alreadyAllocated_sgemm=false;

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
  if (transa=='T' || transa=='t') 
    op_t1 = rocblas_operation_transpose;
  if (transb=='T' || transb=='t') 
    op_t2 = rocblas_operation_transpose;

  rocblas_handle handle = get_rocblas_handle();

  if (!alreadyAllocated_sgemm){
    hipHostMalloc(&Aarray_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Barray_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Carray_sgemm,batchCount*sizeof(float*));
  }
  alreadyAllocated_sgemm=true;

  hipMalloc(&d_Aarray_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Barray_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Carray_sgemm,batchCount*sizeof(float*));

  for (int i = 0; i < batchCount; i++) {
    Aarray_sgemm[i] = (float*) &(A[i*lda*tda]);
    Barray_sgemm[i] = (float*) &(B[i*ldb*tdb]);
    Carray_sgemm[i] = (float*) &(C[i*ldc*tdc]);
  }
  hipMemcpy(d_Aarray_sgemm,Aarray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_sgemm,Barray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_sgemm,Carray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);

  ROCBLAS_CHECK(rocblas_sgemm_batched(handle, op_t1, op_t2, m, n, k,
                                      &alpha,
                                      (const float**) d_Aarray_sgemm, lda,
                                      (const float**) d_Barray_sgemm, ldb,
                                      &beta,
                                      (float**) d_Carray_sgemm, ldc,
                                      batchCount));

  hipDeviceSynchronize();
  
  hipFree(d_Aarray_sgemm);
  hipFree(d_Barray_sgemm);
  hipFree(d_Carray_sgemm);
}

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