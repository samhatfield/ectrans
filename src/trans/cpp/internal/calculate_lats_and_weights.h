// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef CALCULATE_LATS_AND_WEIGHTS_H
#define CALCULATE_LATS_AND_WEIGHTS_H

#include <span>

// Calculate Gaussian latitudes (actually their sine) and weights for Gauss-Legendre quadrature
void calculate_lats_and_weights(std::span<double> mu, std::span<double> weights);

#endif // CALCULATE_LATS_AND_WEIGHTS_H
