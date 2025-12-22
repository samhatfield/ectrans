// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef CONSTANTS_H
#define CONSTANTS_H

// Class for storing physical constants
class Constants {
  private:
    Constants() = default;
    ~Constants() = default;
    inline static Constants* instance = nullptr;

    // Members
    double earth_radius;

  public:
    static Constants& get_instance() {
      if (!instance) {
        instance = new Constants();
      }
      return *instance;
    }

  Constants(const Constants&) = delete("Copy construction disabled");
  Constants& operator=(const Constants&) = delete("Copy assignment disabled");

  static void destroy() {
    delete instance;
    instance = nullptr;
  }

  [[nodiscard]]
  double get_earth_radius() const noexcept {
    return earth_radius;
  }

  void set_earth_radius(double _earth_radius) noexcept {
    earth_radius = _earth_radius;
  }
};

#endif // CONSTANTS_H
