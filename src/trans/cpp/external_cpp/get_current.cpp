// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

// -------------------------------------------------------------------------------------------------
// Fortran binding
// -------------------------------------------------------------------------------------------------

extern "C" {
  // The optional output arguments are passed as pointers which are null when the corresponding
  // argument was not present in the Fortran call.
  void get_current(int* kresol, bool* klam) {
  }
}

// -------------------------------------------------------------------------------------------------
