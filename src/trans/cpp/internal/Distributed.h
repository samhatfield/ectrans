// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include <vector>
#include <span>

// Class for storing resolution-specific distribution parameters
class Distributed {
  private:
    std::vector<int> my_ms; // all of the zonal wavenumbers my W set is responsible for

  public:
    Distributed(int truncation) noexcept {
      my_ms.reserve(truncation + 1);
      for (int m = 0; m < truncation + 1; ++m) {
        my_ms.push_back(m);
      }
    }

    [[nodiscard]] int get_nump() const noexcept { return my_ms.size(); }
    [[nodiscard]] std::span<const int> get_my_ms() const noexcept { return my_ms; }
};

#endif // DISTRIBUTED_H
