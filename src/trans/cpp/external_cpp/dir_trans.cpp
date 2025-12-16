// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dir_trans_sp(
    float* pspvor, float* pspdiv,
    float* pspscalar, float* pspsc3a, float* pspsc3b, float* pspsc2,
    bool ldlatlon, int kproma,
    const int* kvsetuv, const int* kvsetsc,
    int kresol,
    const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
    const float* pgp,
    const float* pgpuv, const float* pgp3a, const float* pgp3b, const float * pgp2) {

    dir_trans(
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      ldlatlon, kproma,
      kvsetuv, kvsetsc,
      kresol,
      kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }

  void dir_trans_dp(
    double* pspvor, double* pspdiv,
    double* pspscalar, double* pspsc3a, double* pspsc3b, double* pspsc2,
    bool ldlatlon, int kproma,
    const int* kvsetuv, const int* kvsetsc,
    int kresol,
    const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
    const double* pgp,
    const double* pgpuv, const double* pgp3a, const double* pgp3b, const double * pgp2) {

    dir_trans(
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      ldlatlon, kproma,
      kvsetuv, kvsetsc,
      kresol,
      kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }
}

// -------------------------------------------------------------------------------------------------
