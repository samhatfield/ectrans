// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef FIND_POLY_ROOTS_H
#define FIND_POLY_ROOTS_H

#include <span>

// 
void find_poly_roots(
    std::span<double> legpol_four, double* lat, double* weight, int num_latitudes, int* iter,
    double* mod
);

#endif // FIND_POLY_ROOTS_H
