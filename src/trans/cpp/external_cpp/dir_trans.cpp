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

#include "abor1.h"

const int nprtrv = 1;
const int mysetv = 1;

// Struct for storing metadata about arguments to DIR_TRANS
struct args_info {
  int ispvor_shape[2];
  int ispdiv_shape[2];
  int ispscalar_shape[2];
  int ispsc3a_shape[3];
  int ispsc3b_shape[3];
  int ispsc2_shape[2];
  int ivsetuv_shape[1];
  int ivsetsc_shape[1];
  int ivsetsc3a_shape[1];
  int ivsetsc3b_shape[1];
  int ivsetsc2_shape[1];
  int igp_shape[3];
  int igpuv_shape[4];
  int igp3a_shape[4];
  int igp3b_shape[4];
  int igp2_shape;
};

template <typename Real> void dir_trans(
  args_info* args,
  Real* pspvor, Real* pspdiv,
  Real* pspscalar, Real* pspsc3a, Real* pspsc3b, Real* pspsc2,
  bool ldlatlon, int kproma,
  const int* kvsetuv, const int* kvsetsc,
  int kresol,
  const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
  const Real* pgp,
  const Real* pgpuv, const Real* pgp3a, const Real* pgp3b, const Real * pgp2) {

  int if_uv = 0;
  int if_uv_g = 0;
  int if_scalars = 0;
  int if_scalars_g = 0;

  if (kvsetuv) {
    // Get total number of UV fields
    if_uv_g = args->ivsetuv_shape[0];

    // Determine which ones are mine
    for (int j = 0; j < if_uv_g; ++j) {
      if (kvsetuv[j] > nprtrv || kvsetuv[j] < 1) {
        std::cerr << "dir_trans: kvsetuv(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetuv too long or contains values outside range");
      }
      if (kvsetuv[j] == mysetv) {
        if_uv += 1;
      }
    }
  } else {
    // No V-set decomposition -> all fields resident on a single V set
    if_uv = args->ispvor_shape[0];
    if_uv_g = if_uv;
  }

  if (kvsetsc) {
    // Get total number of scalar fields
    if_scalars_g = args->ivsetsc_shape[0];

    // Determine which ones are mine
    for (int j = 0; j < if_scalars_g; ++j) {
      if (kvsetsc[j] > nprtrv || kvsetsc[j] < 1) {
        std::cerr << "dir_trans: kvsetsc(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetsc too long or contains values outside range");
      }
      if (kvsetsc[j] == mysetv) {
        if_scalars += 1;
      }
    }
  } else {
    // No V-set decomposition -> all fields resident on a single V set
    if_scalars = args->ispscalar_shape[0];
    if_scalars_g = if_scalars;
  }

  std::cout << "This rank has " << if_uv << " of " << if_uv_g << " UV fields" << std::endl;
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dir_trans_sp(
    args_info* args,
    float* pspvor, float* pspdiv,
    float* pspscalar, float* pspsc3a, float* pspsc3b, float* pspsc2,
    bool ldlatlon, int kproma,
    const int* kvsetuv, const int* kvsetsc,
    int kresol,
    const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
    const float* pgp,
    const float* pgpuv, const float* pgp3a, const float* pgp3b, const float * pgp2) {

    dir_trans(
      args,
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
    args_info* args,
    double* pspvor, double* pspdiv,
    double* pspscalar, double* pspsc3a, double* pspsc3b, double* pspsc2,
    bool ldlatlon, int kproma,
    const int* kvsetuv, const int* kvsetsc,
    int kresol,
    const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
    const double* pgp,
    const double* pgpuv, const double* pgp3a, const double* pgp3b, const double * pgp2) {

    dir_trans(
      args,
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
