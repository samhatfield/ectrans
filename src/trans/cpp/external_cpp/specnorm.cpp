// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to SPECNORM
struct args_info {
  int ipnorm_shape[1];
  int ipspec_shape[2];
  int ikvset_shape[1];
  int ipmet_shape[1];
};
}  // namespace

// The optional scalar arguments are passed as pointers which are null when the corresponding
// argument was not present in the Fortran call.
template <typename Real> void specnorm(
    args_info* args,
    const int* kmaster, const int* kresol,
    Real* pnorm, const Real* pspec, const int* kvset, const Real* pmet) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void specnorm_sp(
      args_info* args,
      const int* kmaster, const int* kresol,
      float* pnorm, const float* pspec, const int* kvset, const float* pmet) {

    specnorm(args, kmaster, kresol, pnorm, pspec, kvset, pmet);
  }

  void specnorm_dp(
      args_info* args,
      const int* kmaster, const int* kresol,
      double* pnorm, const double* pspec, const int* kvset, const double* pmet) {

    specnorm(args, kmaster, kresol, pnorm, pspec, kvset, pmet);
  }
}

// -------------------------------------------------------------------------------------------------
