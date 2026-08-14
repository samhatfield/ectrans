// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef CALCULATE_LATS_AND_WEIGHTS_H
#define CALCULATE_LATS_AND_WEIGHTS_H

#include <vector>

// Gaussian latitudes (actually their sine) and weights for Gauss-Legendre quadrature
struct GaussianQuadrature {
  std::vector<double> mu; // Sine of latitudes
  std::vector<double> weights; // Gaussian weights
};

// Calculate Gaussian latitudes (actually their sine) and weights for Gauss-Legendre quadrature
[[nodiscard]] GaussianQuadrature calculate_lats_and_weights(int num_latitudes);

#endif // CALCULATE_LATS_AND_WEIGHTS_H
