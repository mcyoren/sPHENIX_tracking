#pragma once

#include <cstddef>
#include <vector>

class TpcPolyHelixFitter
{
 public:
  struct Point
  {
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double adc {1.0};
  };

  struct Seed
  {
    bool valid {false};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double px {0.0};
    double py {0.0};
    double pz {0.0};
  };

  struct Options
  {
    // Enforce a deterministic TPC-like ordering before angle unwrapping.
    bool sort_inside_out {true};

    // ADC weights are normalized to the median ADC and clipped.
    bool use_adc_weights {true};
    double adc_weight_power {0.5};
    double min_adc_weight {0.5};
    double max_adc_weight {2.0};

    // Robust iterative fit. Extreme points are removed one at a time.
    bool reject_outliers {true};
    int max_outlier_iterations {3};
    double huber_threshold {2.5};
    double outlier_threshold {5.0};

    // Floors prevent tiny numerical residuals from producing fake outliers.
    double min_sigma_xy {0.02};
    double min_sigma_z {0.05};

    std::size_t min_points {3};
  };

  struct FitResult
  {
    bool ok {false};
    double d0 {0.0};
    double z0 {0.0};
    double phi0 {0.0};
    double theta {0.0};
    double curvature {0.0};
    double xc {0.0};
    double yc {0.0};
    double radius {0.0};
    double dzds {0.0};
    double chi2_xy {0.0};
    double chi2_z {0.0};
    int ndof_xy {0};
    int ndof_z {0};
    int npoints_used {0};
    int npoints_rejected {0};
  };

  static bool fit(const std::vector<Point>& points, FitResult& fit);
  static bool fit(const std::vector<Point>& points, const Seed& seed, FitResult& fit);
  static bool fit(const std::vector<Point>& points, const Seed& seed,
                  FitResult& fit, const Options& options);
};
