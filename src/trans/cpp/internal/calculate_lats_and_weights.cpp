// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include "calculate_lats_and_weights.h"
#include <vector>
#include <cmath>

void calculate_lats_and_weights(std::span<double> mu, std::span<double> weights) {
  // Calculate Gaussian latitudes (actually their sine) and weights for Gauss-Legendre quadrature
  int num_latitudes = static_cast<int>(mu.size());
  int n1 = num_latitudes + 1;
  std::vector<double> legpol_four(n1 * n1);
  double intermediate;
  int odd;

  // Fourier coefficients of series expansion for the ordinary Legendre polynomials 
  legpol_four[0] = 2.0;
  for (int i = 1; i <= num_latitudes; ++i) {
    intermediate = legpol_four[0];
    for (int j = 1; j <= i; ++j) {
      intermediate *= sqrt(1.0 - 0.25 / (double)(j * j));
    }

    odd = i % 2;
    legpol_four[i + n1 * i] = intermediate;
    for (int j = 2; j <= i - odd; j += 2) {
      legpol_four[i + n1 * (i - j)] = legpol_four[i + n1 * (i - j + 2)] * 
        (double)((j - 1) * (2 * i - j + 2)) / (double)(j * (2 * i - j + 1));
    }
  }

  std::vector<double> fn(num_latitudes / 2 + 1);

  double pi = 2.0 * std::asin(1.0);
  odd = num_latitudes % 2;
  int k = odd;
  for (int i = odd; i <= num_latitudes; i+=2) {
    fn[k] = legpol_four[num_latitudes + n1 * i];
    k++;
  }

  std::vector<double> lats(num_latitudes), reg(num_latitudes), li(num_latitudes);

  for (int i = 0; i < num_latitudes / 2; ++i) {
    double z = (double)(4 * (i + 1) - 1) * pi / (double)(4 * num_latitudes + 2);
    lats[i] = z + 1.0 / (tan(z) * (double)(8 * num_latitudes * num_latitudes));
    reg[i] = cos(z);
    li[i] = cos(lats[i]);
  }

  // Refine
  int iter;
  double mod;
  for (int i = num_latitudes / 2; i >= 1; --i) {
    find_poly_roots(fn, &lats[i], &weights[i], num_latitudes, &iter, &mod);
  }

  for (int i = 0; i < num_latitudes / 2; ++i) {
    mu[i] = cos(lats[i]);
  }
}
