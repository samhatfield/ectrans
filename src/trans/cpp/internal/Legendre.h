// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef LEGENDRE_H
#define LEGENDRE_H

template <typename Real> class Legendre {
  private:
    Real* mu; // Sine of latitudes
    Real* weights; // Gaussian weights
    Real* mu_sq_r; // Cosine squared of latitude
    Real* mu_sq_r_sqrt_r; // 1 over cosine of latitude
    Real* epsi; // Epsilon values for the Legendre transform

  public:
    Legendre() {
    }

    [[nodiscard]] int get_trunc() const noexcept { return trunc; }
    [[nodiscard]] int get_num_lats() const noexcept { return num_lats; }
    [[nodiscard]] int get_max_num_lons() const noexcept { return max_num_lons; }
    [[nodiscard]] int get_num_spec_els_glob() const noexcept { return num_spec_els_glob; }
    [[nodiscard]] int get_num_spec_els_glob2() const noexcept { return num_spec_els_glob2; }
    [[nodiscard]] int get_num_nh_lats() const noexcept { return num_nh_lats; }
    [[nodiscard]] int* get_num_lons_per_lat() const noexcept { return num_lons_per_lat; }
};

#endif  // LEGENDRE_H