// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include "calculate_lats_and_weights.h"
#include "find_poly_roots.h"

void find_poly_roots(
    std::span<double> legpol_four, double* lat, double* weight, int num_latitudes, int* iter, 
    double* mod
) {

  int iter_max = 20;
  double x = *lat;
  int odd = num_latitudes % 2;

  for (int i = 0; i <= iter_max; ++i) {
    *iter = i;

    dlk = 0.0;
    if (odd == 0) dlk = 0.5 * legpol_four[0];
    dlxn = 0.0;
    dlldn = 0.0;
    k = 1;

    for (int n = 2 - odd; n <= num_latitudes; n += 2) {
      dlk += fn[k] * cos((double)n * dlx);
      dlldn -= fn[k] * (double)n * sin((double)n * x);
      k++;
    }
    dlmod = -dlk / dlldn;
    dlxn = x + dlmod;
    *mod = dlmod;

    x = 
  }

  *lat = 


}