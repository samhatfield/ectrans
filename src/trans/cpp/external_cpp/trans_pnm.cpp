// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to TRANS_PNM
struct args_info {
  int iprpnm_shape[2];
};
}  // namespace

// The optional scalar/logical arguments are passed as pointers which are null when the
// corresponding argument was not present in the Fortran call.
template <typename Real> void trans_pnm(
    args_info* args,
    const int* kresol, int km, const bool* ktranspose, const bool* kcheap,
    Real* prpnm) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void trans_pnm_sp(
      args_info* args,
      const int* kresol, int km, const bool* ktranspose, const bool* kcheap,
      float* prpnm) {

    trans_pnm(args, kresol, km, ktranspose, kcheap, prpnm);
  }

  void trans_pnm_dp(
      args_info* args,
      const int* kresol, int km, const bool* ktranspose, const bool* kcheap,
      double* prpnm) {

    trans_pnm(args, kresol, km, ktranspose, kcheap, prpnm);
  }
}

// -------------------------------------------------------------------------------------------------
