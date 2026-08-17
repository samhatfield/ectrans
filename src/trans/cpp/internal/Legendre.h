// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <vector>
#include <iostream>
#include "calculate_lats_and_weights.h"
#include "LegendrePolynomialGroup.h"

template <typename Real> class Legendre {
  private:
    GaussianQuadrature quadrature; // Sine of latitudes and Gaussian weights
    std::vector<double> mu_sq_r; // Cosine squared of latitude
    std::vector<double> mu_sq_r_sqrt_r; // 1 over cosine of latitude
    std::vector<double> epsi; // Epsilon values for the Legendre transform
    LegendrePolynomialGroup<Real> legendre_polynomial_group; // Group of Legendre polynomials

  public:
    Legendre(
      int num_latitudes, int truncation, std::span<const int> my_ms, int nprtrv,
      std::span<const int> num_latitudes_m
    ) : quadrature(calculate_lats_and_weights(num_latitudes)),
        mu_sq_r(num_latitudes),
        mu_sq_r_sqrt_r(num_latitudes),
        epsi(num_latitudes),
        legendre_polynomial_group(my_ms, nprtrv, truncation, num_latitudes, num_latitudes_m) {
    }

    [[nodiscard]] double* get_mu() noexcept { return quadrature.mu.data(); }
    [[nodiscard]] double* get_weights() noexcept { return quadrature.weights.data(); }
    [[nodiscard]] double* get_mu_sq_r() noexcept { return mu_sq_r.data(); }
    [[nodiscard]] double* get_mu_sq_r_sqrt_r() noexcept { return mu_sq_r_sqrt_r.data(); }
    [[nodiscard]] double* get_epsi() noexcept { return epsi.data(); }
    [[nodiscard]] LegendrePolynomialGroup<Real>& get_legendre_polynomial_group() noexcept {
      return legendre_polynomial_group;
    }
};

#endif  // LEGENDRE_H
