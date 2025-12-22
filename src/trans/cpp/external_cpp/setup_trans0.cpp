// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <iostream>

#include "abor1.h"
#include "General.h"
#include "Constants.h"
#include "Distributed.h"

// -------------------------------------------------------------------------------------------------
// Fortran binding
// -------------------------------------------------------------------------------------------------

extern "C" {
  void setup_trans0(
    int kout, int kerr, int kprintlev, int kmax_resol, int kpromatr, int kprgpns, int kprgpew,
    int kprtrw, int kmpoff, int ksync_trans, int ktrans_sync_level, int keq_regions, double prad,
    int kalloperm, int kopt_memory_tr, int* k_regions_ns, int* k_regions_ew, int* k_regions) {

    // Default values
    int nout = 6;
    int nerr = 0;
    int nprintlev = 0;
    int nmax_resol = 1;
    int nprgpns = 1;
    int nprgpew = 1;
    int nprtrw = 1;
    int n_regions_ns = 1;
    int n_regions_ew = 1;
    int npromatr = 0;
    bool lmpoff = false;
    bool lsync_trans = false;
    int ntrans_sync_level = 0;
    int leq_regions = false;
    double ra = 6371229.0;
    bool lalloperm = false;
    int nstack_memory_tr = 0;

    if (kout >= 0) nout = kout;
    if (kerr >= 0) nerr = kerr;
    if (kprintlev >= 0) nprintlev = kprintlev;

    if (nprintlev > 0) std::cout << "Entering routine setup_trans0" << std::endl;

    if (kmax_resol >= 0) nmax_resol = kmax_resol;
    if (kpromatr >= 0) {
      if (kpromatr % 2 != 0) ABOR1("setup_trans0: kpromatr must be a multiple of 2");
      npromatr = kpromatr;
    }
    if (kprgpns >= 0) nprgpns = kprgpns;
    if (kprgpew >= 0) nprgpew = kprgpew;
    if (kprtrw >= 0) nprtrw = kprtrw;
    if (kmpoff >= 0) lmpoff = kmpoff;
    if (ksync_trans >= 0) lsync_trans = ksync_trans;
    if (ktrans_sync_level >= 0) ntrans_sync_level = ktrans_sync_level;
    if (keq_regions >= 0) leq_regions = keq_regions;
    if (kalloperm >= 0) lalloperm = kalloperm;
    if (kopt_memory_tr >= 0) nstack_memory_tr = kopt_memory_tr;

    General::get_instance().init(
      nout, nerr, nprintlev, nmax_resol, npromatr, lalloperm, lmpoff, lsync_trans,
      ntrans_sync_level, kopt_memory_tr
    );

    // Set Earth radius
    double default_earth_radius = 6371229.0;
    if (prad > 0.0) {
      Constants::get_instance().set_earth_radius(prad);
    } else {
      Constants::get_instance().set_earth_radius(default_earth_radius);
    }

    // Initialise resolution-agnostic parallelisation parameters
    Distributed::get_instance().init(nprgpns, nprgpew, nprtrw, leq_regions);

    if (*k_regions_ns >= 0) *k_regions_ns = Distributed::get_instance().get_n_regions_ns();
    if (*k_regions_ew >= 0) *k_regions_ew = Distributed::get_instance().get_n_regions_ew();
    if (k_regions) k_regions = Distributed::get_instance().get_n_regions();
  }
}

// -------------------------------------------------------------------------------------------------
