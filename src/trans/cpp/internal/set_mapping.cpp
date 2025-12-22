// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include "set_mapping.h"
#include "abor1.h"

// Map task (1-indexed) to NS/EW and W/V set
void task2sets(
  int task, bool leq_regions, int nprgpew, int nproc, int nprtrv, int* n_regions, int n_regions_ns,
  int* my_region_ns, int* my_region_ew, int* mysetw, int* mysetv
) {
  if (task <= 0 || task > nproc) {
    ABOR1("task2set: invalid argument");
  } else {
    if (leq_regions) {
      *my_region_ns = 1;
      int _task = task;
      for (int i = 0; i < n_regions_ns; i++) {
        if (_task > n_regions[i] ) {
          _task -= n_regions[i];
          *my_region_ns += 1;
          continue;
        }
        *my_region_ew = _task;
        break;
      }
    } else {
      *my_region_ew = (task - 1) % nprgpew + 1;
      *my_region_ns = (task - 1) / nprgpew + 1;
    }

    *mysetv = (task - 1) % nprtrv + 1;
    *mysetw = (task - 1) / nprtrv + 1;
  }
}
