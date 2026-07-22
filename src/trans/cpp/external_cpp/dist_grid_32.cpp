// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to DIST_GRID_32
struct args_info {
  int ipgpg_shape[2];
  int ikfrom_shape[1];
  int igp_shape[3];
};
}  // namespace

// PGPG and PGP are always single precision. The optional scalar arguments are passed as pointers
// which are null when the corresponding argument was not present in the Fortran call.
void dist_grid_32_impl(
    args_info* args,
    const int* kproma, int kfdistg, const int* kresol,
    const float* pgpg, const int* kfrom, float* pgp) {
}

// -------------------------------------------------------------------------------------------------
// Fortran binding
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dist_grid_32(
      args_info* args,
      const int* kproma, int kfdistg, const int* kresol,
      const float* pgpg, const int* kfrom, float* pgp) {

    dist_grid_32_impl(args, kproma, kfdistg, kresol, pgpg, kfrom, pgp);
  }
}

// -------------------------------------------------------------------------------------------------
