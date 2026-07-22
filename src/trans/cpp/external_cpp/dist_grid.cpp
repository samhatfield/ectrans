// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to DIST_GRID
struct args_info {
  int ipgpg_shape[2];
  int ikfrom_shape[1];
  int igp_shape[3];
  int iksort_shape[1];
};
}  // namespace

// The optional scalar arguments are passed as pointers which are null when the corresponding
// argument was not present in the Fortran call.
template <typename Real> void dist_grid(
    args_info* args,
    const int* kproma, int kfdistg, const int* kresol,
    const Real* pgpg, const int* kfrom, Real* pgp, const int* ksort) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dist_grid_sp(
      args_info* args,
      const int* kproma, int kfdistg, const int* kresol,
      const float* pgpg, const int* kfrom, float* pgp, const int* ksort) {

    dist_grid(args, kproma, kfdistg, kresol, pgpg, kfrom, pgp, ksort);
  }

  void dist_grid_dp(
      args_info* args,
      const int* kproma, int kfdistg, const int* kresol,
      const double* pgpg, const int* kfrom, double* pgp, const int* ksort) {

    dist_grid(args, kproma, kfdistg, kresol, pgpg, kfrom, pgp, ksort);
  }
}

// -------------------------------------------------------------------------------------------------
