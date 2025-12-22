// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef GENERAL_H
#define GENERAL_H

// Class for storing general configuration parameters
class General {
  private:
    General() = default;
    ~General() = default;
    inline static General* instance = nullptr;

    // Members
    int outstream;
    int errstream;
    int print_level;
    int max_resolutions;
    int npromatr;
    bool permanent_allocations;
    bool mp_off;
    bool synchronization;
    int synchronization_level;
    int transpose_mem_strategy;

  public:
    static General& get_instance() {
      if (!instance) {
        instance = new General();
      }
      return *instance;
    }

  General(const General&) = delete("Copy construction disabled");
  General& operator=(const General&) = delete("Copy assignment disabled");

  static void destroy() {
    delete instance;
    instance = nullptr;
  }

  [[nodiscard]] int get_outstream() const noexcept { return outstream; }
  [[nodiscard]] int get_errstream() const noexcept { return errstream; }
  [[nodiscard]] int get_print_level() const noexcept { return print_level; }
  [[nodiscard]] int get_max_resolutions() const noexcept { return max_resolutions; }
  [[nodiscard]] int get_npromatr() const noexcept { return npromatr; }
  [[nodiscard]] bool use_permanent_allocations() const noexcept { return permanent_allocations; }
  [[nodiscard]] bool mp_is_off() const noexcept { return mp_off; }
  [[nodiscard]] bool use_synchronization() const noexcept { return synchronization; }
  [[nodiscard]] int get_synchronization_level() const noexcept { return synchronization_level; }
  [[nodiscard]] int get_transpose_mem_strategy() const noexcept { return transpose_mem_strategy; }

  void init(
    int _outstream, int _errstream, int _print_level, int _max_resolutions, int _npromatr,
    bool _permanent_allocations, bool _mp_off, bool _synchronization, int _synchronization_level,
    int _transpose_mem_strategy
  ) noexcept {
    outstream = _outstream;
    errstream = _errstream;
    print_level = _print_level;
    max_resolutions = _max_resolutions;
    npromatr = _npromatr;
    permanent_allocations = _permanent_allocations;
    mp_off = _mp_off;
    synchronization = _synchronization;
    synchronization_level = _synchronization_level;
    transpose_mem_strategy = _transpose_mem_strategy;
  }
};

#endif // GENERAL_H
