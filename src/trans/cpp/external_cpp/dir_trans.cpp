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

//  struct args_info{
//    int* ispvor_shape, ispdiv_shape, ispscalar_shape, ispsc3a_shape, ispsc3b_shape, ispsc2_shape;
//    int* ivsetuv_shape, ivsetsc_shape, ivsetsc3a_shape, ivsetsc3b_shape, ivsetsc2_shape;
//    int* igp_shape, igpuv_shape, igp3a_shape, igp3b_shape, igp2_shape;
//  };

struct args_info{
  int ispvor_shape[2];
  int ivsetuv_shape[1];
};

template <typename Real> void dir_trans(
  args_info args,
  Real* pspvor, Real* pspdiv,
  Real* pspscalar, Real* pspsc3a, Real* pspsc3b, Real* pspsc2,
  bool ldlatlon, int kproma,
  const int* kvsetuv, const int* kvsetsc,
  int kresol,
  const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
  const Real* pgp,
  const Real* pgpuv, const Real* pgp3a, const Real* pgp3b, const Real * pgp2) {

  if (std::is_same<Real, float>::value) {
    std::cout << "Inside dir_trans_sp" << std::endl;
  } else {
    std::cout << "Inside dir_trans_dp" << std::endl;
  }

  if (kvsetuv) {
    std::cout << " kvsetuv is present" << std::endl;
    // std::cout << args.ispvor_shape[0] << std::endl;
    // // std::cout << args.ispvor_shape[1] << std::endl;
    // std::cout << args.ivsetuv_shape[0] << std::endl;
    std::cout << " done" << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dir_trans_sp(
    args_info args,
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
    args_info args,
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
