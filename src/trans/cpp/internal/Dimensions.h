// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include <algorithm>
#include <vector>
#include <span>
#include "abor1.h"

class Dimensions {
  private:
    int trunc; // Spectral truncation
    int num_lats; // Number of latitudes pole-to-pole
    int max_num_lons; // Maximum number of longitudes on any latitude
    int num_spec_els_glob; // Number of complex spectral coefficients
    int num_spec_els_glob2; // Number of complex spectral coefficients * 2
    int num_nh_lats; // Number of latitudes pole-to-equator
    int* num_lons_per_lat; // Number of longitudes on each latitude
    bool reduced_grid; // True if the grid is reduced, false if it is regular
    std::vector<int> num_latitudes_m; // Number of latitudes for each zonal wavenumber

  public:
    Dimensions(int _trunc, int _num_lats, int* _max_num_lons, int* _num_lons_per_lat) {
      trunc = _trunc;
      num_lats = _num_lats;

      if (num_lats <= 0 || num_lats % 2 != 0) ABOR1("Dimensions: num_lats not positive and even");

      if (_max_num_lons) {
        max_num_lons = *_max_num_lons;
      } else {
        max_num_lons = 2 * num_lats;
      }

      num_lons_per_lat = new int[num_lats];
      if (_num_lons_per_lat) {
        reduced_grid = _num_lons_per_lat[0] != _num_lons_per_lat[1];
        max_num_lons = 0;
        for (int i = 0; i < num_lats; ++i) {
          if (_num_lons_per_lat[i] <= 0) {
            ABOR1("Dimensions.init: must have positive definite number of longitudes per latitude");
          }
          max_num_lons = std::max(max_num_lons, _num_lons_per_lat[i]);
          num_lons_per_lat[i] = _num_lons_per_lat[i];
        }
      } else {
        for (int i = 0; i < num_lats; ++i) num_lons_per_lat[i] = max_num_lons;
      }

      num_spec_els_glob = (trunc + 1) * (trunc + 2) / 2;
      num_spec_els_glob2 = 2 * num_spec_els_glob;
      num_nh_lats = num_lats / 2;

      // Compute number of latitudes considered for the Legendre transform at each zonal wavenumber
      num_latitudes_m.reserve(trunc + 1);
      for (int m = 0; m <= trunc; ++m) {
        // TODO: Make this actually dependent on latitude
        num_latitudes_m.push_back(num_lats);
      }
    }

    [[nodiscard]] int get_trunc() const noexcept { return trunc; }
    [[nodiscard]] int get_num_lats() const noexcept { return num_lats; }
    [[nodiscard]] int get_max_num_lons() const noexcept { return max_num_lons; }
    [[nodiscard]] int get_num_spec_els_glob() const noexcept { return num_spec_els_glob; }
    [[nodiscard]] int get_num_spec_els_glob2() const noexcept { return num_spec_els_glob2; }
    [[nodiscard]] int get_num_nh_lats() const noexcept { return num_nh_lats; }
    [[nodiscard]] int* get_num_lons_per_lat() const noexcept { return num_lons_per_lat; }
    [[nodiscard]] std::span<const int> get_num_latitudes_m() const noexcept { return num_latitudes_m; }
};

#endif // DIMENSIONS_H
