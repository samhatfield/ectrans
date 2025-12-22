// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include "abor1.h"
#include "General.h"
#include "set_mapping.h"

// Class for storing parallel configuration parameters
class Distributed {
  private:
    Distributed() = default;
    ~Distributed() = default;
    inline static Distributed* instance = nullptr;

    // Members
    int nproc;
    int nprgpns;
    int nprgpew;
    int nprtrw;
    int nprtrv;
    int nprtrns;
    bool leq_regions;
    int myproc;
    int* n_regions;
    int n_regions_ns;
    int n_regions_ew;
    int my_region_ns;
    int my_region_ew;
    int mysetw;
    int mysetv;
    int* nprcids;

  public:
    static Distributed& get_instance() {
      if (!instance) {
        instance = new Distributed();
      }
      return *instance;
    }

  Distributed(const Distributed&) = delete("Copy construction disabled");
  Distributed& operator=(const Distributed&) = delete("Copy assignment disabled");

  static void destroy() {
    delete instance;
    instance = nullptr;
  }

  [[nodiscard]] int get_nproc() const noexcept { return nproc; }
  [[nodiscard]] int get_nprgpns() const noexcept { return nprgpns; }
  [[nodiscard]] int get_nprgpew() const noexcept { return nprgpew; }
  [[nodiscard]] int get_nprtrw() const noexcept { return nprtrw; }
  [[nodiscard]] int get_nprtrv() const noexcept { return nprtrv; }
  [[nodiscard]] int get_nprtrns() const noexcept { return nprtrns; }
  [[nodiscard]] bool use_equal_regions() const noexcept { return leq_regions; }
  [[nodiscard]] int get_myproc() const noexcept { return myproc; }
  [[nodiscard]] int* get_n_regions() const noexcept { return n_regions; }
  [[nodiscard]] int get_n_regions_ns() const noexcept { return n_regions_ns; }
  [[nodiscard]] int get_n_regions_ew() const noexcept { return n_regions_ew; }
  [[nodiscard]] int get_my_region_ns() const noexcept { return my_region_ns; }
  [[nodiscard]] int get_my_region_ew() const noexcept { return my_region_ew; }
  [[nodiscard]] int get_mysetw() const noexcept { return mysetw; }
  [[nodiscard]] int get_mysetv() const noexcept { return mysetv; }
  [[nodiscard]] int* get_nprcids() const noexcept { return nprcids; }

  void init(int _nprgpns, int _nprgpew, int _nprtrw, int _leq_regions) noexcept {
    nprgpns = _nprgpns;
    nprgpew = _nprgpew;
    nproc = nprgpns * nprgpew;
    nprtrw = _nprtrw;
    nprtrns = nprtrw;

    if (nproc % nprtrw != 0 || nprtrw > nproc) {
      ABOR1("Distributed.init: nproc inconsistent with nprtrw");
    }

    nprtrv = nproc / nprtrw;

    if (General::get_instance().get_print_level() > 0) {
      std::cout << "nproc = " << nproc << std::endl;
      std::cout << "nprgpns = " << nprgpns << std::endl;
      std::cout << "nprgpew = " << nprgpew << std::endl;
      std::cout << "nprtrw = " << nprtrw << std::endl;
      std::cout << "nprtrv = " << nprtrv << std::endl;
    }

    if (nproc > 1) {
      // Not supported yet
      ABOR1("Not implemented");
    } else {
      myproc = 1;
    }

    if (myproc > nproc) {
      ABOR1("Distributed.init: myproc > nproc");
    }

    leq_regions = _leq_regions;
    if (leq_regions) {
      ABOR1("Not implemented");
    } else {
      n_regions_ns = nprgpns;
      n_regions = new int[n_regions_ns];
      for (int i = 0; i < n_regions_ns; ++i) n_regions[i] = nprgpew;
      n_regions_ew = nprgpew;
    }

    task2sets(
      myproc, leq_regions, nprgpew, nproc, nprtrv, n_regions, n_regions_ns,
      &my_region_ns, &my_region_ew, &mysetw, &mysetv
    );

    if (General::get_instance().get_print_level() > 0) {
      std::cout << "myproc = " << myproc << std::endl;
      std::cout << "my_region_ns = " << my_region_ns << std::endl;
      std::cout << "my_region_ew = " << my_region_ew << std::endl;
      std::cout << "mysetw = " << mysetw << std::endl;
      std::cout << "mysetv = " << mysetv << std::endl;
    }

    nprcids = new int[nproc];
    for (int i = 0; i < nproc; ++i) nprcids[i] = i;

    // TODO: create communicators for MPI groups
  }
};

#endif // DISTRIBUTED_H
