// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef DIMENSIONS_H
#define DIMENSIONS_H

class Dimensions {
  private:
    int trunc; // Spectral truncation
    int num_lats; // Number of latitudes pole-to-pole
    int max_num_lons; // Maximum number of longitudes on any latitude
    int num_spec_els_glob; // Number of complex spectral coefficients
    int num_spec_els_glob2; // Number of complex spectral coefficients * 2
    int num_nh_lats; // Number of latitudes pole-to-equator

  public:
    Dimensions(int _trunc, int _num_lats, int _max_num_lons) {
      trunc = _trunc;
      num_lats = _num_lats;

      if (num_lats <= 0 || num_lats % 2 != 0) ABOR1("Dimensions: num_lats not positive and even");

      if (_max_num_lons > 0) {
        max_num_lons = _max_num_lons;
      } else {
        max_num_lons = 2 * num_lats;
      }

      num_spec_els_glob = (trunc + 1) * (trunc + 2) / 2;
      num_spec_els_glob2 = 2 * num_spec_els_glob2;
      num_nh_lats = (num_lats + 1) / 2;
    }

    [[nodiscard]] int get_trunc() const noexcept { return trunc; }
    [[nodiscard]] int get_num_lats() const noexcept { return num_lats; }
};

#endif // DIMENSIONS_H
