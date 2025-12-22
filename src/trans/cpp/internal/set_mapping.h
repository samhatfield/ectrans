// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef SET_MAPPING_H
#define SET_MAPPING_H

// Map task (1-indexed) to NS/EW and W/V set
void task2sets(
  int task, bool leq_regions, int nprgpew, int nproc, int nprtrv, int* n_regions, int n_regions_ns,
  int* my_region_ns, int* my_region_ew, int* mysetw, int* mysetv
);

#endif // SET_MAPPING_H
