// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to GPNORM_TRANS
struct args_info {
  int igp_shape[3];
  int ipave_shape[1];
  int ipmin_shape[1];
  int ipmax_shape[1];
};
}  // namespace

// KRESOL is optional and passed as a pointer which is null when the argument was not present in the
// Fortran call.
template <typename Real> void gpnorm_trans(
    args_info* args,
    int kfields, int kproma, bool kave_only, const int* kresol,
    const Real* pgp, Real* pave, Real* pmin, Real* pmax) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void gpnorm_trans_sp(
      args_info* args,
      int kfields, int kproma, bool kave_only, const int* kresol,
      const float* pgp, float* pave, float* pmin, float* pmax) {

    gpnorm_trans(args, kfields, kproma, kave_only, kresol, pgp, pave, pmin, pmax);
  }

  void gpnorm_trans_dp(
      args_info* args,
      int kfields, int kproma, bool kave_only, const int* kresol,
      const double* pgp, double* pave, double* pmin, double* pmax) {

    gpnorm_trans(args, kfields, kproma, kave_only, kresol, pgp, pave, pmin, pmax);
  }
}

// -------------------------------------------------------------------------------------------------
