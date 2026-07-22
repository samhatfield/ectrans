// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to DIST_SPEC
struct args_info {
  int ipspecg_shape[2];
  int ikfrom_shape[1];
  int ikvset_shape[1];
  int ipspec_shape[2];
  int iksort_shape[1];
};
}  // namespace

// The optional scalar/logical arguments are passed as pointers which are null when the
// corresponding argument was not present in the Fortran call.
template <typename Real> void dist_spec(
    args_info* args,
    const bool* kdim1_is_fld, int kfdistg, const int* kresol, const int* ksmax,
    const Real* pspecg, const int* kfrom, const int* kvset, Real* pspec, const int* ksort) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dist_spec_sp(
      args_info* args,
      const bool* kdim1_is_fld, int kfdistg, const int* kresol, const int* ksmax,
      const float* pspecg, const int* kfrom, const int* kvset, float* pspec, const int* ksort) {

    dist_spec(args, kdim1_is_fld, kfdistg, kresol, ksmax, pspecg, kfrom, kvset, pspec, ksort);
  }

  void dist_spec_dp(
      args_info* args,
      const bool* kdim1_is_fld, int kfdistg, const int* kresol, const int* ksmax,
      const double* pspecg, const int* kfrom, const int* kvset, double* pspec, const int* ksort) {

    dist_spec(args, kdim1_is_fld, kfdistg, kresol, ksmax, pspecg, kfrom, kvset, pspec, ksort);
  }
}

// -------------------------------------------------------------------------------------------------
