// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef RESOLUTION_DEPENDENT_H
#define RESOLUTION_DEPENDENT_H

#include <vector>

#include "abor1.h"
#include "Dimensions.h"
#include "Distributed.h"
#include "Legendre.h"
#include "Parallel.h"

// Class for managing resolution-dependent structures
template <typename Real> class ResolutionDependent {
  private:
    ResolutionDependent() = default;
    ~ResolutionDependent() = default;
    inline static ResolutionDependent* instance = nullptr;

    // Members
    std::vector<Dimensions> dimension_list;
    std::vector<Distributed> distributed_list;
    std::vector<Legendre<Real>> legendre_list;

  public:
    static ResolutionDependent& get_instance() {
      if (!instance) {
        instance = new ResolutionDependent();
      }
      return *instance;
    }

  ResolutionDependent(const ResolutionDependent&) = delete;
  ResolutionDependent& operator=(const ResolutionDependent&) = delete;

  static void destroy() {
    delete instance;
    instance = nullptr;
  }

  [[nodiscard]]
  Dimensions& get_dimension_list(int resol) {
    if (resol <= dimension_list.size()) {
      return dimension_list[resol - 1];
    } else {
      ABOR1("ResolutionDependent.get_dimensions: resol requested does not exist");
    }
  }

  [[nodiscard]]
  Distributed& get_distributed_list(int resol) {
    if (resol <= distributed_list.size()) {
      return distributed_list[resol - 1];
    } else {
      ABOR1("ResolutionDependent.get_distributed: resol requested does not exist");
    }
  }

  [[nodiscard]]
  Legendre<Real>& get_legendre_list(int resol) {
    if (resol <= legendre_list.size()) {
      return legendre_list[resol - 1];
    } else {
      ABOR1("ResolutionDependent.get_legendre: resol requested does not exist");
    }
  }

  int init_resol(
    int truncation, int num_latitudes, int* max_lons_per_lat, int* num_lons_per_lat
  ) noexcept {
    // Dimensions (truncation, number of latitudes etc.)
    dimension_list.push_back(Dimensions(
      truncation, num_latitudes, max_lons_per_lat, num_lons_per_lat
    ));

    // Resolution-dependent distribution parameters (number of zonal wavenumbers per W set etc.)
    distributed_list.push_back(Distributed(truncation));

    // Legendre transform-related arrays (Gaussian weights, Legendre polynomials etc.)
    legendre_list.push_back(Legendre<Real>(
      num_latitudes, truncation, distributed_list.back().get_my_ms(),
      Parallel::get_instance().get_nprtrv(), dimension_list.back().get_num_latitudes_m()
    ));

    // Return index integer handle (like KRESOL in Fortran, so it starts from 1)
    return dimension_list.size();
  }

  void delete_resol(int resol) noexcept {
    dimension_list.erase(dimension_list.begin() + resol - 1);
    distributed_list.erase(distributed_list.begin() + resol - 1);
    legendre_list.erase(legendre_list.begin() + resol - 1);
  }
};

#endif // RESOLUTION_DEPENDENT_H
