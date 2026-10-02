#include "DCGemCurrent.h"

R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)

#include <iostream>
#include <iomanip>

void test_dc_gem_current(int run = 81558,
                         ULong64_t bco = 0,
                         const char* cdbGlobalTag = "newcdbtag")
{
  if (bco == 0)
  {
    std::cerr << "Provide a representative DC BCO, e.g. from dc_windows_summary_run<run>.csv" << std::endl;
    return;
  }

  const auto r = dcgem::load_r1_currents(run, bco, cdbGlobalTag);

  constexpr ULong64_t GTM_MASK = (1ULL << 40U) - 1ULL;

  std::cout << std::setprecision(16)
            << "available=" << r.available << '\n'
            << "dc_target_bco_40bit=" << (r.target_bco & GTM_MASK) << '\n'
            << "aligned_target_bco_absolute=" << r.aligned_target_bco << '\n'
            << "first_conditions_bco_absolute=" << r.first_conditions_bco << '\n'
            << "last_conditions_bco_absolute=" << r.last_conditions_bco << '\n'
            << "first_conditions_bco_40bit=" << (r.first_conditions_bco & GTM_MASK) << '\n'
            << "last_conditions_bco_40bit=" << (r.last_conditions_bco & GTM_MASK) << '\n'
            << "nearest_conditions_bco_absolute=" << r.nearest_bco << '\n'
            << "LoadSR1=" << r.load_sr1 << '\n'
            << "LoadNR1=" << r.load_nr1 << '\n'
            << "AverageSR1=" << r.average_sr1 << '\n'
            << "AverageNR1=" << r.average_nr1 << '\n'
            << "CDB URL=" << r.cdb_url << std::endl;
}
