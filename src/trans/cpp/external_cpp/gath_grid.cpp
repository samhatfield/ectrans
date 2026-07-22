// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to GATH_GRID
struct args_info {
  int ipgpg_shape[2];
  int ikto_shape[1];
  int igp_shape[3];
};
}  // namespace

// The optional scalar arguments are passed as pointers which are null when the corresponding
// argument was not present in the Fortran call.
template <typename Real> void gath_grid(
    args_info* args,
    const int* kproma, int kfgathg, const int* kresol,
    Real* pgpg, const int* kto, const Real* pgp) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void gath_grid_sp(
      args_info* args,
      const int* kproma, int kfgathg, const int* kresol,
      float* pgpg, const int* kto, const float* pgp) {

    gath_grid(args, kproma, kfgathg, kresol, pgpg, kto, pgp);
  }

  void gath_grid_dp(
      args_info* args,
      const int* kproma, int kfgathg, const int* kresol,
      double* pgpg, const int* kto, const double* pgp) {

    gath_grid(args, kproma, kfgathg, kresol, pgpg, kto, pgp);
  }
}

// -------------------------------------------------------------------------------------------------
