// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include "calculate_lats_and_weights.h"
#include "find_poly_roots.h"
#include <limits>
#include <cmath>

void find_poly_roots(
    std::span<double> legpol_four, double* lat, double* weight, int num_latitudes, int* iter, 
    double* mod
) {

  static constexpr int iter_max = 20;
  double x = *lat;
  int odd = num_latitudes % 2;

  static constexpr double eps = std::numeric_limits<double>::epsilon();

  for (int i = 1; i <= iter_max + 1; ++i) {
    double dlldn = 0.0;
    int k = 1;
    if (abs(*mod) <= eps * 1000.0) {
      // Last pass
      for (int j = 2 - odd; j <= num_latitudes; j += 2) {
        // Normalised derivative
        dlldn -= legpol_four[k] * (double)j * sin((double)j * (*lat));
        k++;
      }
      *weight = (double)(2 * num_latitudes + 1) / (dlldn * dlldn);
      return;
    }

    double dlk = 0.0;
    if (odd == 0) dlk = 0.5 * legpol_four[0];

    for (int j = 2 - odd; j <= num_latitudes; j += 2) {
      // Normalised ordinary Legendre polynomial == \overbar{P_n}^0
      dlk += legpol_four[k] * cos((double)j * (*lat));
      // Normalised derivative == d/d\theta(\overbar{P_n}^0)
      dlldn -= legpol_four[k] * (double)j * sin((double)j * (*lat));
      k++;
    }

    // Newton method
    *mod = -dlk / dlldn;
    *lat += *mod;
  }
}
