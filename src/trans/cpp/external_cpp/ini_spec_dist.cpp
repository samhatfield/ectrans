// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to INI_SPEC_DIST
struct args_info {
  int ikasm0_shape[1];
  int ikprocm_shape[1];
  int ikumpp_shape[1];
  int ikpossp_shape[1];
  int ikmyms_shape[1];
  int ikptrms_shape[1];
  int ikallms_shape[1];
};
}  // namespace

// -------------------------------------------------------------------------------------------------
// Fortran binding
// -------------------------------------------------------------------------------------------------

extern "C" {
  // The optional output arguments are passed as pointers which are null when the corresponding
  // argument was not present in the Fortran call.
  void ini_spec_dist(
      args_info* args,
      int ksmax, int ktmax, int kprtrw, int kmysetw,
      int* kspolegl, int* kspec, int* kspec2, int* kspec2mx,
      int* kasm0, int* kprocm, int* kumpp, int* kpossp, int* kmyms, int* kptrms, int* kallms) {
  }
}

// -------------------------------------------------------------------------------------------------
