// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Struct for storing metadata (shapes) about the array arguments to INV_TRANS
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

// The optional scalar/logical arguments are passed as pointers which are null when the
// corresponding argument was not present in the Fortran call.
template <typename Real> void inv_trans(
    args_info* args,
    const bool* kscders, const bool* kvorgp, const bool* kdivgp, const bool* kuvder,
    const bool* kdlatlon, const int* kproma, const int* kresol,
    const Real* pspvor, const Real* pspdiv,
    const Real* pspscalar, const Real* pspsc3a, const Real* pspsc3b, const Real* pspsc2,
    const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
    const int* kvsetsc2,
    Real* pgp,
    Real* pgpuv, Real* pgp3a, Real* pgp3b, Real* pgp2) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void inv_trans_sp(
      args_info* args,
      const bool* kscders, const bool* kvorgp, const bool* kdivgp, const bool* kuvder,
      const bool* kdlatlon, const int* kproma, const int* kresol,
      const float* pspvor, const float* pspdiv,
      const float* pspscalar, const float* pspsc3a, const float* pspsc3b, const float* pspsc2,
      const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
      const int* kvsetsc2,
      float* pgp,
      float* pgpuv, float* pgp3a, float* pgp3b, float* pgp2) {

    inv_trans(
      args,
      kscders, kvorgp, kdivgp, kuvder, kdlatlon, kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc, kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }

  void inv_trans_dp(
      args_info* args,
      const bool* kscders, const bool* kvorgp, const bool* kdivgp, const bool* kuvder,
      const bool* kdlatlon, const int* kproma, const int* kresol,
      const double* pspvor, const double* pspdiv,
      const double* pspscalar, const double* pspsc3a, const double* pspsc3b, const double* pspsc2,
      const int* kvsetuv, const int* kvsetsc, const int* kvsetsc3a, const int* kvsetsc3b,
      const int* kvsetsc2,
      double* pgp,
      double* pgpuv, double* pgp3a, double* pgp3b, double* pgp2) {

    inv_trans(
      args,
      kscders, kvorgp, kdivgp, kuvder, kdlatlon, kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc, kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }
}

// -------------------------------------------------------------------------------------------------
