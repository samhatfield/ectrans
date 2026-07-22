// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to VORDIV_TO_UV
struct args_info {
  int ipspvor_shape[2];
  int ipspdiv_shape[2];
  int ipspu_shape[2];
  int ipspv_shape[2];
  int ikvsetuv_shape[1];
};
}  // namespace

// KVSETUV is optional and passed as a pointer which is null when the argument was not present in
// the Fortran call.
template <typename Real> void vordiv_to_uv(
    args_info* args,
    int ksmax,
    const Real* pspvor, const Real* pspdiv, Real* pspu, Real* pspv, const int* kvsetuv) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void vordiv_to_uv_sp(
      args_info* args,
      int ksmax,
      const float* pspvor, const float* pspdiv, float* pspu, float* pspv, const int* kvsetuv) {

    vordiv_to_uv(args, ksmax, pspvor, pspdiv, pspu, pspv, kvsetuv);
  }

  void vordiv_to_uv_dp(
      args_info* args,
      int ksmax,
      const double* pspvor, const double* pspdiv, double* pspu, double* pspv, const int* kvsetuv) {

    vordiv_to_uv(args, ksmax, pspvor, pspdiv, pspu, pspv, kvsetuv);
  }
}

// -------------------------------------------------------------------------------------------------
