// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to GPNORM_TRANSTL
struct args_info {
  int igp_shape[3];
  int ipave_shape[1];
};
}  // namespace

// KRESOL is optional and passed as a pointer which is null when the argument was not present in the
// Fortran call.
template <typename Real> void gpnorm_transtl(
    args_info* args,
    int kfields, int kproma, const int* kresol,
    const Real* pgp, Real* pave) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void gpnorm_transtl_sp(
      args_info* args,
      int kfields, int kproma, const int* kresol,
      const float* pgp, float* pave) {

    gpnorm_transtl(args, kfields, kproma, kresol, pgp, pave);
  }

  void gpnorm_transtl_dp(
      args_info* args,
      int kfields, int kproma, const int* kresol,
      const double* pgp, double* pave) {

    gpnorm_transtl(args, kfields, kproma, kresol, pgp, pave);
  }
}

// -------------------------------------------------------------------------------------------------
