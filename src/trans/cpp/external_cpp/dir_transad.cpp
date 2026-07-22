// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Struct for storing metadata (shapes) about the array arguments to DIR_TRANSAD
struct args_info {
  int ispvor_shape[2];
  int ispdiv_shape[2];
  int ispscalar_shape[2];
  int ispsc3a_shape[3];
  int ispsc3b_shape[3];
  int ispsc2_shape[2];
  int ivsetuv_shape[1];
  int ivsetsc_shape[1];
  int ivsetsc3a_shape[1];
  int ivsetsc3b_shape[1];
  int ivsetsc2_shape[1];
  int igp_shape[3];
  int igpuv_shape[4];
  int igp3a_shape[4];
  int igp3b_shape[4];
  int igp2_shape[3];
};
}  // namespace

// The optional scalar arguments are passed as pointers which are null when the corresponding
// argument was not present in the Fortran call.
template <typename Real> void dir_transad(
    args_info* args,
    const int* kproma, const int* kresol,
    Real* pspvor, Real* pspdiv,
    Real* pspscalar, Real* pspsc3a, Real* pspsc3b, Real* pspsc2,
    const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
    const int* kvsetsc2,
    Real* pgp,
    Real* pgpuv, Real* pgp3a, Real* pgp3b, Real* pgp2) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dir_transad_sp(
      args_info* args,
      const int* kproma, const int* kresol,
      float* pspvor, float* pspdiv,
      float* pspscalar, float* pspsc3a, float* pspsc3b, float* pspsc2,
      const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
      const int* kvsetsc2,
      float* pgp,
      float* pgpuv, float* pgp3a, float* pgp3b, float* pgp2) {

    dir_transad(
      args,
      kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc, kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }

  void dir_transad_dp(
      args_info* args,
      const int* kproma, const int* kresol,
      double* pspvor, double* pspdiv,
      double* pspscalar, double* pspsc3a, double* pspsc3b, double* pspsc2,
      const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
      const int* kvsetsc2,
      double* pgp,
      double* pgpuv, double* pgp3a, double* pgp3b, double* pgp2) {

    dir_transad(
      args,
      kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc, kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }
}

// -------------------------------------------------------------------------------------------------
