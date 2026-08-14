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
#include "Legendre.h"

// Class for managing resolution-dependent structures
template <typename Real> class ResolutionDependent {
  private:
    ResolutionDependent() = default;
    ~ResolutionDependent() = default;
    inline static ResolutionDependent* instance = nullptr;

    // Members
    std::vector<Dimensions> dimension_list;
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
    dimension_list.push_back(Dimensions(truncation, num_latitudes, max_lons_per_lat, num_lons_per_lat));
    legendre_list.push_back(Legendre<Real>(num_latitudes));

    return dimension_list.size();
  }

  void delete_resol(int resol) noexcept {
    dimension_list.erase(dimension_list.begin() + resol - 1);
    legendre_list.erase(legendre_list.begin() + resol - 1);
  }
};

#endif // RESOLUTION_DEPENDENT_H
