// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include <vector>
#include <span>
#include "Parallel.h"
#include "setup_spectral_distribution.h"

// Class for storing resolution-specific distribution parameters
class Distributed {
  private:
    std::vector<int> moffset; // starting index of each zonal wavenumber in spectral arrays
    int spolegl; // TODO: not sure what this is
    std::vector<int> mtask; // which W-set owns each zonal wavenumber
    std::vector<int> num_ms; // how many zonal wavenumbers each W-set owns
    int num_my_ms; // how many zonal wavenumbers my W-set owns
    int nspec; // how many complex spectral coefficients this task owns
    int nspec2; // nspec * 2
    int nspec2max; // maximum number of complex spectral coefficients on any W-set * 2
    std::vector<int> wsetoffset; // starting index of each W-set in global spectral arrays
    std::vector<int> my_ms; // all of the zonal wavenumbers my W set is responsible for

  public:
    Distributed(int truncation) noexcept {
      std::vector<int> my_ms_all_ms(truncation + 1); // extra large version of my_ms
      num_ms.reserve(Parallel::get_instance().get_nprtrw());
      int mysetw = Parallel::get_instance().get_mysetw();
      setup_spectral_distribution(
        truncation, Parallel::get_instance().get_nprtrw(), mysetw,
        num_ms, my_ms_all_ms
      );
      num_my_ms = num_ms[mysetw - 1];
      my_ms.assign(my_ms_all_ms.begin(), my_ms_all_ms.begin() + num_my_ms);
    }

    [[nodiscard]] int get_nump() const noexcept { return my_ms.size(); }
    [[nodiscard]] std::span<const int> get_my_ms() const noexcept { return my_ms; }
};

#endif // DISTRIBUTED_H
