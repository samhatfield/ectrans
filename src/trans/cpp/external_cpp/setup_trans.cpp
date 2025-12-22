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

template <typename Real> void setup_trans(int ksmax, int kdgl, int kdlon, int* kresol) {
  int resolution_handle = ResolutionDependent::get_instance().init_resol(ksmax, kdgl, kdlon);

  if (General::get_instance().get_print_level() > 0) {
    std::cout << "Defined resolution " << resolution_handle << std::endl;
  }

  if (*kresol == 0) {
    *kresol = resolution_handle;
  }
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void setup_trans_sp(int ksmax, int kdgl, int kdlon, int* kresol) {
    setup_trans<float>(ksmax, kdgl, kdlon, kresol);
  }

  void setup_trans_dp(int ksmax, int kdgl, int kdlon, int* kresol) {
    setup_trans<double>(ksmax, kdgl, kdlon, kresol);
  }
}

// -------------------------------------------------------------------------------------------------
