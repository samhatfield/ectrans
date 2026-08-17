// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <cstring>

#include "Constants.h"
#include "Distributed.h"
#include "General.h"
#include "ResolutionDependent.h"
#include "kokkos_lifecycle.h"

// -------------------------------------------------------------------------------------------------
// Fortran binding
// -------------------------------------------------------------------------------------------------

extern "C" {
  // cdmode is a null-terminated C string, or null when the optional argument was not present in the
  // Fortran call.
  void trans_end(const char* cdmode) {
    // As in the CPU backend, the mode defaults to 'FINAL' when the argument is absent. 'INTER' asks
    // for a partial teardown that keeps the library usable, so it must not tear down the runtime.
    const char* mode = cdmode ? cdmode : "FINAL";

    if (std::strcmp(mode, "FINAL") != 0) return;

    // Order matters. Every Kokkos::View owned by ecTrans lives inside one of these singletons, and
    // Kokkos raises an error if an allocation is freed after Kokkos::finalize, so the singletons must
    // be destroyed first. Both precisions are destroyed because this routine is compiled once into
    // ectrans_cpp_common and shared by the single- and double-precision libraries; destroying a
    // singleton that was never instantiated is a no-op.
    ResolutionDependent<float>::destroy();
    ResolutionDependent<double>::destroy();
    General::destroy();
    Constants::destroy();

    finalise_kokkos();
  }
}

// -------------------------------------------------------------------------------------------------
