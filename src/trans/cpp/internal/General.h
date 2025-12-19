// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

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
    bool synchronisation;
    int synchronisation_level;
    int transposition_memory_strategy;

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

  [[nodiscard]]
  double get_outstream() const noexcept {
    return outstream;
  }

  [[nodiscard]]
  double get_errstream() const noexcept {
    return errstream;
  }

  void init(
    int _outstream, int _errstream, int _print_level, int _max_resolutions, int _npromatr,
    bool _permanent_allocations, bool _mp_off, bool _synchronisation, int _synchronisation_level,
    int _transposition_memory_strategy
  ) noexcept {
    outstream = _outstream;
  }
};
