// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef LEGENDRE_POLYNOMIAL_GROUP_H
#define LEGENDRE_POLYNOMIAL_GROUP_H

#include <Kokkos_Core.hpp>

template <typename Real> class LegendrePolynomialGroup {
  private:
    std::vector<Kokkos::View<Real**>> anti;
    std::vector<Kokkos::View<Real**>> symm;

  public:
    LegendrePolynomialGroup(
        std::span<const int> my_ms, int nprtrv, int truncation, int num_latitudes,
        std::span<const int> num_latitudes_m
    )
      : anti(my_ms.size()), symm(my_ms.size()) {
      for (int mloc = 0; mloc < my_ms.size(); ++mloc) {
        int prtrv = std::min(nprtrv, static_cast<int>(my_ms.size()) - mloc);

        for (int this_mloc = mloc; this_mloc < mloc + prtrv; ++this_mloc) {
          int m = my_ms[this_mloc];
          int ila = (truncation - m + 2) / 2;
          int idglu = std::min(num_latitudes, num_latitudes_m[m]);

          anti[mloc] = Kokkos::View<Real**>("anti", idglu, ila);
          std::cout << "Allocated anti[" << mloc << "] with dimensions (" << idglu << ", " << ila << ")" << std::endl;
        }
      }
    }
};

#endif // LEGENDRE_POLYNOMIAL_GROUP_H
