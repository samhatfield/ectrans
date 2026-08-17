// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef KOKKOS_LIFECYCLE_H
#define KOKKOS_LIFECYCLE_H

#include <iostream>

#include <Kokkos_Core.hpp>

#include "abor1.h"

// ecTrans is a library, so it is not necessarily the owner of the Kokkos runtime. The calling
// application (the IFS, Atlas, or a standalone driver) may well have called Kokkos::initialize
// itself, in which case it also owns Kokkos allocations that outlive trans_end, and finalising the
// runtime from underneath it would abort the run. We therefore record whether ecTrans was the one
// that started Kokkos, and only finalise what we started.
//
// This is an inline variable rather than a member of a singleton because the flag must be readable
// before any ecTrans object exists, and shared by every translation unit that includes this header.
// None of this is thread-safe: setup_trans0 and trans_end are documented as collective, single-
// threaded calls made once per run.
inline bool ectrans_initialised_kokkos = false;

// Bring up the Kokkos runtime, unless it is already running. Must be called before anything that
// allocates a Kokkos::View, which in practice means at the top of setup_trans0. Calling it more than
// once is harmless. When verbose is true, Kokkos prints its configuration (backends, architecture,
// thread count) during initialisation.
inline void initialise_kokkos(bool verbose) {
  if (Kokkos::is_initialized()) {
    if (verbose) {
      std::cout << "initialise_kokkos: Kokkos is already initialised, ecTrans will use the existing"
                << " runtime and will not finalise it in trans_end" << std::endl;
    }
    return;
  }

  // Kokkos cannot be restarted once it has been finalised, so there is nothing sensible to do here
  // but fail loudly. This is reached if the calling application finalises Kokkos before calling
  // setup_trans0, or if a trans_end('FINAL') is followed by another setup_trans0.
  if (Kokkos::is_finalized()) {
    ABOR1("initialise_kokkos: Kokkos has already been finalised and cannot be reinitialised");
  }

  // The settings are populated programmatically rather than from argc/argv, which the Fortran
  // interface does not give us. Anything not set explicitly here is still picked up from the
  // KOKKOS_* environment variables, so the usual command-line tuning remains available.
  Kokkos::InitializationSettings settings;
  settings.set_print_configuration(verbose);

  Kokkos::initialize(settings);
  ectrans_initialised_kokkos = true;

  if (verbose) {
    std::cout << "initialise_kokkos: initialised Kokkos with default execution space "
              << Kokkos::DefaultExecutionSpace::name() << std::endl;
  }
}

// Shut the Kokkos runtime down, but only if ecTrans started it. Every Kokkos::View owned by ecTrans
// must already have been destroyed by the time this is called, because Kokkos raises an error when
// an allocation is freed after Kokkos::finalize. See the ordering comment in trans_end.
inline void finalise_kokkos() {
  // Either the calling application owns the runtime, or setup_trans0 was never reached. Either way
  // finalisation is not ours to do.
  if (!ectrans_initialised_kokkos) return;

  // Defend against the calling application having finalised Kokkos in the meantime. Kokkos::finalize
  // is not idempotent and aborts on a second call.
  if (!Kokkos::is_initialized()) {
    ectrans_initialised_kokkos = false;
    return;
  }

  Kokkos::finalize();
  ectrans_initialised_kokkos = false;
}

#endif // KOKKOS_LIFECYCLE_H
