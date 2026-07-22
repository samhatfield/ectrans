// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdlib.h>

namespace {
// Shapes of the array arguments to TRANS_INQ
struct args_info {
  int ikgptotl_shape[2];
  int ikmyms_shape[1];
  int ikasm0_shape[1];
  int ikumpp_shape[1];
  int ikpossp_shape[1];
  int ikptrms_shape[1];
  int ikallms_shape[1];
  int ikdim0g_shape[1];
  int ikfrstlat_shape[1];
  int iklstlat_shape[1];
  int ikptrlat_shape[1];
  int ikptrfrstlat_shape[1];
  int ikptrlstlat_shape[1];
  int iksta_shape[2];
  int ikonl_shape[2];
  int ikultpp_shape[1];
  int ikptrls_shape[1];
  int iknmeng_shape[1];
  int ildsplitlat_shape[1];
  int iplapin_shape[1];
  int iknvalue_shape[1];
  int ipmu_shape[1];
  int ipgw_shape[1];
  int iprpnm_shape[2];
  int ikpms_shape[1];
  int ikdglu_shape[1];
};
}  // namespace

// KRESOL is an input; all the other arguments are outputs. Each optional argument is passed as a
// pointer which is null when the corresponding argument was not present in the Fortran call. PMU is
// always double precision.
template <typename Real> void trans_inq(
    args_info* args,
    const int* kresol,
    int* kspec, int* kspec2, int* kspec2g, int* kspec2mx, int* knump,
    int* kgptot, int* kgptotg, int* kgptotmx, int* kfrstloff, int* kptrfloff,
    int* kprtrw, int* kmysetw, int* kmysetv, int* kmy_region_ns, int* kmy_region_ew,
    int* klei3, int* kspolegl, int* ksmax, int* kdef_resol, bool* klam,
    int* kgptotl, int* kmyms, int* kasm0, int* kumpp, int* kpossp, int* kptrms, int* kallms,
    int* kdim0g, int* kfrstlat, int* klstlat, int* kptrlat, int* kptrfrstlat, int* kptrlstlat,
    int* ksta, int* konl, int* kultpp, int* kptrls, int* knmeng, bool* ksplitlat,
    Real* plapin, int* knvalue, double* pmu, Real* pgw, Real* prpnm, int* kpms, int* kdglu) {
}

// -------------------------------------------------------------------------------------------------
// Fortran bindings
// -------------------------------------------------------------------------------------------------

extern "C" {
  void trans_inq_sp(
      args_info* args,
      const int* kresol,
      int* kspec, int* kspec2, int* kspec2g, int* kspec2mx, int* knump,
      int* kgptot, int* kgptotg, int* kgptotmx, int* kfrstloff, int* kptrfloff,
      int* kprtrw, int* kmysetw, int* kmysetv, int* kmy_region_ns, int* kmy_region_ew,
      int* klei3, int* kspolegl, int* ksmax, int* kdef_resol, bool* klam,
      int* kgptotl, int* kmyms, int* kasm0, int* kumpp, int* kpossp, int* kptrms, int* kallms,
      int* kdim0g, int* kfrstlat, int* klstlat, int* kptrlat, int* kptrfrstlat, int* kptrlstlat,
      int* ksta, int* konl, int* kultpp, int* kptrls, int* knmeng, bool* ksplitlat,
      float* plapin, int* knvalue, double* pmu, float* pgw, float* prpnm, int* kpms, int* kdglu) {

    trans_inq(
      args, kresol, kspec, kspec2, kspec2g, kspec2mx, knump, kgptot, kgptotg, kgptotmx, kfrstloff,
      kptrfloff, kprtrw, kmysetw, kmysetv, kmy_region_ns, kmy_region_ew, klei3, kspolegl, ksmax,
      kdef_resol, klam, kgptotl, kmyms, kasm0, kumpp, kpossp, kptrms, kallms, kdim0g, kfrstlat,
      klstlat, kptrlat, kptrfrstlat, kptrlstlat, ksta, konl, kultpp, kptrls, knmeng, ksplitlat,
      plapin, knvalue, pmu, pgw, prpnm, kpms, kdglu
    );
  }

  void trans_inq_dp(
      args_info* args,
      const int* kresol,
      int* kspec, int* kspec2, int* kspec2g, int* kspec2mx, int* knump,
      int* kgptot, int* kgptotg, int* kgptotmx, int* kfrstloff, int* kptrfloff,
      int* kprtrw, int* kmysetw, int* kmysetv, int* kmy_region_ns, int* kmy_region_ew,
      int* klei3, int* kspolegl, int* ksmax, int* kdef_resol, bool* klam,
      int* kgptotl, int* kmyms, int* kasm0, int* kumpp, int* kpossp, int* kptrms, int* kallms,
      int* kdim0g, int* kfrstlat, int* klstlat, int* kptrlat, int* kptrfrstlat, int* kptrlstlat,
      int* ksta, int* konl, int* kultpp, int* kptrls, int* knmeng, bool* ksplitlat,
      double* plapin, int* knvalue, double* pmu, double* pgw, double* prpnm, int* kpms,
      int* kdglu) {

    trans_inq(
      args, kresol, kspec, kspec2, kspec2g, kspec2mx, knump, kgptot, kgptotg, kgptotmx, kfrstloff,
      kptrfloff, kprtrw, kmysetw, kmysetv, kmy_region_ns, kmy_region_ew, klei3, kspolegl, ksmax,
      kdef_resol, klam, kgptotl, kmyms, kasm0, kumpp, kpossp, kptrms, kallms, kdim0g, kfrstlat,
      klstlat, kptrlat, kptrfrstlat, kptrlstlat, ksta, konl, kultpp, kptrls, knmeng, ksplitlat,
      plapin, knvalue, pmu, pgw, prpnm, kpms, kdglu
    );
  }
}

// -------------------------------------------------------------------------------------------------
