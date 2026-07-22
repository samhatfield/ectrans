// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <iostream>

#include "General.h"
#include "ResolutionDependent.h"

// kdlon and kresol are optional and passed as pointers which are null when the corresponding
// argument was not present in the Fortran call.
template <typename Real> void setup_trans(int ksmax, int kdgl, int* kdlon, int* kloen, int* kresol) {
  int resolution_handle = ResolutionDependent::get_instance().init_resol(ksmax, kdgl, kdlon, kloen);

  if (General::get_instance().get_print_level() > 0) {
    std::cout << "Defined resolution " << resolution_handle << std::endl;
  }

  if (kresol) {
    *kresol = resolution_handle;
  }
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void setup_trans_sp(int ksmax, int kdgl, int* kdlon, int* kloen, int* kresol) {
    setup_trans<float>(ksmax, kdgl, kdlon, kloen, kresol);
  }

  void setup_trans_dp(int ksmax, int kdgl, int* kdlon, int* kloen, int* kresol) {
    setup_trans<double>(ksmax, kdgl, kdlon, kloen, kresol);
  }
}

// -------------------------------------------------------------------------------------------------
