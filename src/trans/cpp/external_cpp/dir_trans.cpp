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
int ngptot = 100;
int nproma;
bool latlon;
int ngpblks;

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
  int igp2_shape[3];
};

template <typename Real> void dir_trans(
    args_info* args,
    int kdlatlon, int kproma, int kresol,
    Real* pspvor, Real* pspdiv,
    Real* pspscalar, Real* pspsc3a, Real* pspsc3b, Real* pspsc2,
    const int* kvsetuv, const int* kvsetsc,
    const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
    const Real* pgp,
    const Real* pgpuv, const Real* pgp3a, const Real* pgp3b, const Real * pgp2) {

  int if_uv = 0;
  int if_uv_g = 0;
  int if_scalars = 0;
  int if_scalars_g = 0;
  int nf_sc2 = 0;
  int nf_sc3a = 0;
  int nf_sc3b = 0;
  int if_sc2_g = 0;
  int if_sc3a_g = 0;
  int if_sc3b_g = 0;
  nproma = ngptot;
  latlon = false;

  if (kvsetuv) {
    // Get total number of UV fields
    if_uv_g = args->ivsetuv_shape[0];

    // Determine which ones are mine
    for (int j = 0; j < if_uv_g; ++j) {
      if (kvsetuv[j] > nprtrv || kvsetuv[j] < 1) {
        std::cerr << "dir_trans: kvsetuv(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetuv too long or contains values outside range");
      }
      if (kvsetuv[j] == mysetv) if_uv += 1;
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
      if (kvsetsc[j] == mysetv) if_scalars += 1;
    }
  } else {
    // No V-set decomposition -> all fields resident on a single V set
    if_scalars = args->ispscalar_shape[0];
    if_scalars_g = if_scalars;
  }

  if (kvsetsc2) {
    // Get total number of 2D scalar fields
    if (!pspsc2) {
      ABOR1("dir_trans: kvsetsc2 present but not pspsc2");
    }
    if_sc2_g = args->ivsetsc2_shape[0];
    if_scalars_g += if_sc2_g;
    for (int j = 0; j < args->ivsetsc2_shape[0]; ++j) {
      if (kvsetsc2[j] > nprtrv || kvsetsc2[j] < 1) {
        std::cerr << "dir_trans: kvsetsc2(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetsc2 too long or contains values outside range");
      }
      if (kvsetsc2[j] == mysetv) {
        if_scalars += 1;
        nf_sc2 += 1;
      }
    }
  } else if (pspsc2) {
    if_sc2_g = args->ispsc2_shape[0];
    nf_sc2 = args->ispsc2_shape[0];
    if_scalars += nf_sc2;
    if_scalars_g += if_sc2_g;
  }

  if (kvsetsc3a) {
    // Get total number of 3D scalar fields (3a)
    if (!pspsc3a) ABOR1("dir_trans: kvsetsc3a present but not pspsc3a");
    if_sc3a_g = args->ivsetsc3a_shape[0];
    if_scalars_g += if_sc3a_g * args->ispsc3a_shape[2];
    for (int j = 0; j < args->ivsetsc3a_shape[0]; ++j) {
      if (kvsetsc3a[j] > nprtrv || kvsetsc3a[j] < 1) {
        std::cerr << "dir_trans: kvsetsc3a(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetsc3a too long or contains values outside range");
      }
      if (kvsetsc3a[j] == mysetv) {
        if_scalars += args->ispsc3a_shape[2];
        nf_sc3a += 1;
      }
    }
  } else if (pspsc3a) {
    if_scalars += args->ispsc3a_shape[0] * args->ispsc3a_shape[2];
    if_sc3a_g = args->ispsc3a_shape[0];
    if_scalars_g += if_sc3a_g * args->ispsc3a_shape[2];
    nf_sc3a = args->ispsc3a_shape[0];
  }

  if (kvsetsc3b) {
    // Get total number of 3D scalar fields (3b)
    if (!pspsc3b) ABOR1("dir_trans: kvsetsc3b present but not pspsc3b");
    if_sc3b_g = args->ivsetsc3b_shape[0];
    if_scalars_g += if_sc3b_g * args->ispsc3b_shape[2];
    for (int j = 0; j < args->ivsetsc3b_shape[0]; ++j) {
      if (kvsetsc3b[j] > nprtrv || kvsetsc3b[j] < 1) {
        std::cerr << "dir_trans: kvsetsc3b(" << j << ") > nprtrv or < 1" << std::endl;
        ABOR1("dir_trans: kvsetsc3b too long or contains values outside range");
      }
      if (kvsetsc3b[j] == mysetv) {
        if_scalars += args->ispsc3b_shape[2];
        nf_sc3b += 1;
      }
    }
  } else if (pspsc3b) {
    if_scalars += args->ispsc3b_shape[0] * args->ispsc3b_shape[2];
    if_sc3b_g = args->ispsc3b_shape[0];
    if_scalars_g += if_sc3b_g * args->ispsc3b_shape[2];
    nf_sc3b = args->ispsc3b_shape[0];
  }

  if (kproma >= 0) nproma = kproma;
  if (kdlatlon >= 0) latlon = kdlatlon;

  ngpblks = (ngptot - 1) / nproma + 1;

  int if_fs = 2 * if_uv + if_scalars;
  int if_gp = 2 * if_uv_g + if_scalars_g;

  // Consistency checks
  if (if_uv > 0) {
    if (!pspvor) ABOR1("dir_trans: if_uv > 0 but pspvor missing");
    if (args->ispvor_shape[0] < if_uv) {
      std::cerr << "dir_trans: pspvor too short: " << args->ispvor_shape[0] << " < " << if_uv
                << std::endl;
      ABOR1("dir_trans: pspvor too short");
    }
    if (!pspdiv) ABOR1("dir_trans: pspvor present but pspdiv missing");
    if (args->ispdiv_shape[0] != if_uv) {
      std::cerr << "dir_trans: pspdiv too short: " << args->ispdiv_shape[0] << " < " << if_uv
                << std::endl;
      ABOR1("dir_trans: inconsistent first dimension of pspvor and pspdiv");
    }
  }

  if (if_scalars > 0) {
    if (pspscalar) {
      if (args->ispscalar_shape[0] < if_scalars) {
        std::cerr << "dir_trans: pspscalar too short: " << args->ispscalar_shape[0] << " < "
                  << if_scalars << std::endl;
        ABOR1("dir_trans: pspscalar too short");
      }
      if (pspsc3a || pspsc3b || pspsc2) {
        ABOR1("dir_trans: pspscalar and (pspsc3a or pspsc3b or pspsc2) both present");
      }
    }
  }

  if (nprtrv > 1) {
    if (if_uv > 0 && !kvsetuv) {
      std::cerr << "nprtrv > 1 and if_uv > 0 but kvsetuv missing " << std::endl;
      ABOR1("dir_trans: specify vertical spectral distribution");
    }
    if (pspscalar && !kvsetsc) {
      std::cerr << "nprtrv > 1 and pspscalar present but kvsetsc missing" << std::endl;
      ABOR1("dir_trans: specify vertical spectral distribution");
    }
    if (pspsc2 && !kvsetsc2) {
      std::cerr << "nprtrv > 1 and pspsc2 present but kvsetsc2 missing" << std::endl;
      ABOR1("dir_trans: specify vertical spectral distribution");
    }
    if (pspsc3a && !kvsetsc3a) {
      std::cerr << "nprtrv > 1 and pspsc3a present but kvsetsc3a missing" << std::endl;
      ABOR1("dir_trans: specify vertical spectral distribution");
    }
    if (pspsc3b && !kvsetsc3b) {
      std::cerr << "nprtrv > 1 and pspsc3b present but kvsetsc3b missing" << std::endl;
      ABOR1("dir_trans: specify vertical spectral distribution");
    }
  }

  // Check dimensions of pgp make sense
  if (pgp) {
    if (args->igp_shape[0] < nproma) {
      std::cerr << "dir_trans: first dimension of pgp too small: " << args->igp_shape[0]
                << " < " << nproma << std::endl;
      ABOR1("dir_trans: first dimension of pgp too small");
    }
    if (args->igp_shape[1] < if_gp) {
      std::cerr << "dir_trans: second dimension of pgp too small: " << args->igp_shape[1]
                << " < " << if_gp << std::endl;
      ABOR1("dir_trans: second dimension of pgp too small");
    }
    if (args->igp_shape[2] < ngpblks) {
      std::cerr << "dir_trans: third dimension of pgp too small: " << args->igp_shape[2]
                << " < " << ngpblks << std::endl;
      ABOR1("dir_trans: third dimension of pgp too small");
    }
  }

  if (pgpuv) {
    if (!pspvor) ABOR1("dir_trans: pspvor has to be present when pgpuv is");
    if (args->igpuv_shape[0] < nproma) {
      std::cerr << "dir_trans: first dimension of pgpuv too small: " << args->igpuv_shape[0]
                << " < " << nproma << std::endl;
      ABOR1("dir_trans: first dimension of pgpuv too small");
    }
    if (args->igpuv_shape[1] != if_uv_g) {
      std::cerr << "dir_trans: second dimension of pgpuv inconsistent: " << args->igpuv_shape[1]
                << " != " << if_uv_g << std::endl;
      ABOR1("dir_trans: second dimension of pgpuv inconsistent");
    }
    if (args->igpuv_shape[2] < 2) {
      std::cerr << "dir_trans: third dimension of pgpuv too small: " << args->igpuv_shape[2]
                << " < 2" << std::endl;
      ABOR1("dir_trans: third dimension of pgpuv too small");
    }
    if (args->igpuv_shape[3] < ngpblks) {
      std::cerr << "dir_trans: fourth dimension of pgpuv too small: " << args->igpuv_shape[3]
                << " < " << ngpblks << std::endl;
      ABOR1("dir_trans: fourth dimension of pgpuv too small");
    }
  }

  if (pgp2 && !pspsc2) ABOR1("dir_trans: pspsc2 has to be present when pgp2 is");
  if (if_sc2_g > 0) {
    if (pgp2) {
      if (args->igp2_shape[0] < nproma) {
        std::cerr << "dir_trans: first dimension of pgp2 too small: " << args->igp2_shape[0]
                  << " < " << nproma << std::endl;
        ABOR1("dir_trans: first dimension of pgp2 too small");
      }
      if (args->igp2_shape[1] != if_sc2_g) {
        std::cerr << "dir_trans: second dimension of pgp2 inconsistent: " << args->igp2_shape[1]
                  << " != " << if_sc2_g << std::endl;
        ABOR1("dir_trans: second dimension of pgp2 inconsistent");
      }
      if (args->igp2_shape[2] < ngpblks) {
        std::cerr << "dir_trans: third dimension of pgp2 too small: " << args->igp2_shape[2]
                  << " < " << ngpblks << std::endl;
        ABOR1("dir_trans: third dimension of pgp2 too small");
      }
    }
  }

  if (pgp3a && !pspsc3a) ABOR1("dir_trans: pspsc3a has to be present when pgp3a is");
  if (if_sc3a_g > 0) {
    if (pgp3a) {
      if (args->igp3a_shape[0] < nproma) {
        std::cerr << "dir_trans: first dimension of pgp3a too small: " << args->igp3a_shape[0]
                  << " < " << nproma << std::endl;
        ABOR1("dir_trans: first dimension of pgp3a too small");
      }
      if (args->igp3a_shape[1] != if_sc3a_g) {
        std::cerr << "dir_trans : second dimension of pgp3a inconsistent: " << args->igp3a_shape[1]
                  << " != " << if_sc3a_g << std::endl;
        ABOR1("dir_trans: second dimension of pgp3a inconsistent");
      }
      if (args->igp3a_shape[2] != args->ispsc3a_shape[2]) {
        std::cerr << "dir_trans: third dimension of pgp3a inconsistent: " << args->igp3a_shape[2]
                  << " != " << args->ispsc3a_shape[2] << std::endl;
        ABOR1("dir_trans: third dimension of pgp3a inconsistent");
      }
      if (args->igp3a_shape[3] < ngpblks) {
        std::cerr << "dir_trans: fourth dimension of pgp3a too small: " << args->igp3a_shape[3]
                  << " < " << ngpblks << std::endl;
        ABOR1("dir_trans: fourth dimension of pgp3a too small");
      }
    } else ABOR1("dir_trans: pgp3a missing");
  }

  if (pgp3b && !pspsc3b) ABOR1("dir_trans: pspsc3b has to be present when pgp3b is");
  if (if_sc3b_g > 0) {
    if (pgp3b) {
      if (args->igp3b_shape[0] < nproma) {
        std::cerr << "dir_trans: first dimension of pgp3b too small: " << args->igp3b_shape[0]
                  << " < " << nproma << std::endl;
        ABOR1("dir_trans: first dimension of pgp3b too small");
      }
      if (args->igp3b_shape[1] != if_sc3b_g) {
        std::cerr << "dir_trans : second dimension of pgp3b inconsistent: " << args->igp3b_shape[1]
                  << " != " << if_sc3b_g << std::endl;
        ABOR1("dir_trans: second dimension of pgp3b inconsistent");
      }
      if (args->igp3b_shape[2] != args->ispsc3b_shape[2]) {
        std::cerr << "dir_trans: third dimension of pgp3b inconsistent: " << args->igp3b_shape[2]
                  << " != " << args->ispsc3b_shape[2] << std::endl;
        ABOR1("dir_trans: third dimension of pgp3b inconsistent");
      }
      if (args->igp3b_shape[3] < ngpblks) {
        std::cerr << "dir_trans: fourth dimension of pgp3b too small: " << args->igp3b_shape[3]
                  << " < " << ngpblks << std::endl;
        ABOR1("dir_trans: fourth dimension of pgp3b too small");
      }
    } else ABOR1("dir_trans: pgp3b missing");
  }

  // -----------------------------------------------------------------------------------------------

  // dir_trans_ctl(
  //   if_uv_g, if_scalars_g, if_gp, if_fs, if_uv, if_scalars, pspvor, pspdiv, pspscalar, kvsetuv,
  //   kvsetsc, pgp, pspsc3a, pspsc3b, pspsc2, kvsetsc3a, kvsetsc3b, kvsetsc2, pgpuv, pgp3a, pgp3b,
  //   pgp2
  // )

  // -----------------------------------------------------------------------------------------------
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void dir_trans_sp(
      args_info* args,
      int kdlatlon, int kproma, int kresol,
      float* pspvor, float* pspdiv,
      float* pspscalar, float* pspsc3a, float* pspsc3b, float* pspsc2,
      const int* kvsetuv, const int* kvsetsc,
      const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
      const float* pgp,
      const float* pgpuv, const float* pgp3a, const float* pgp3b, const float * pgp2) {

    dir_trans(
      args,
      kdlatlon, kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc,
      kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }

  void dir_trans_dp(
      args_info* args,
      int kdlatlon, int kproma, int kresol,
      double* pspvor, double* pspdiv,
      double* pspscalar, double* pspsc3a, double* pspsc3b, double* pspsc2,
      const int* kvsetuv, const int* kvsetsc,
      const int* kvsetsc3a, const int* kvsetsc3b, const int* kvsetsc2,
      const double* pgp,
      const double* pgpuv, const double* pgp3a, const double* pgp3b, const double * pgp2) {

    dir_trans(
      args,
      kdlatlon, kproma, kresol,
      pspvor, pspdiv,
      pspscalar, pspsc3a, pspsc3b, pspsc2,
      kvsetuv, kvsetsc,
      kvsetsc3a, kvsetsc3b, kvsetsc2,
      pgp,
      pgpuv, pgp3a, pgp3b, pgp2
    );
  }
}

// -------------------------------------------------------------------------------------------------
