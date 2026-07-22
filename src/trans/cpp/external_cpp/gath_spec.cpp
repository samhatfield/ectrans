// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to GATH_SPEC
struct args_info {
  int ipspecg_shape[2];
  int ikto_shape[1];
  int ikvset_shape[1];
  int ipspec_shape[2];
};
}  // namespace

// The optional scalar/logical arguments are passed as pointers which are null when the
// corresponding argument was not present in the Fortran call.
template <typename Real> void gath_spec(
    args_info* args,
    const bool* kdim1_is_fld, const bool* kza0ip, int kfgathg, const int* kresol, const int* ksmax,
    Real* pspecg, const int* kto, const int* kvset, const Real* pspec) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void gath_spec_sp(
      args_info* args,
      const bool* kdim1_is_fld, const bool* kza0ip, int kfgathg, const int* kresol, const int* ksmax,
      float* pspecg, const int* kto, const int* kvset, const float* pspec) {

    gath_spec(args, kdim1_is_fld, kza0ip, kfgathg, kresol, ksmax, pspecg, kto, kvset, pspec);
  }

  void gath_spec_dp(
      args_info* args,
      const bool* kdim1_is_fld, const bool* kza0ip, int kfgathg, const int* kresol, const int* ksmax,
      double* pspecg, const int* kto, const int* kvset, const double* pspec) {

    gath_spec(args, kdim1_is_fld, kza0ip, kfgathg, kresol, ksmax, pspecg, kto, kvset, pspec);
  }
}

// -------------------------------------------------------------------------------------------------
