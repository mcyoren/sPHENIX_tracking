#include "TpcPolyHelixFitter.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace
{
  using Point = TpcPolyHelixFitter::Point;
  using Seed = TpcPolyHelixFitter::Seed;
  using Options = TpcPolyHelixFitter::Options;
  using FitResult = TpcPolyHelixFitter::FitResult;

  constexpr double two_pi = 2.0 * M_PI;

  struct WorkPoint
  {
    Point point;
    std::size_t original_index {0};
    double base_weight {1.0};
    double robust_weight {1.0};
    bool active {true};
  };

  struct InternalFit
  {
    FitResult result;
    double sign {1.0};
    std::vector<std::size_t> active_indices;
    std::vector<double> s;
    std::vector<double> residual_xy;
    std::vector<double> residual_z;
  };

  double wrap_pi(double phi)
  {
    while (phi > M_PI) phi -= two_pi;
    while (phi <= -M_PI) phi += two_pi;
    return phi;
  }

  double unwrap_near(double phi, const double ref)
  {
    while (phi - ref > M_PI) phi -= two_pi;
    while (phi - ref < -M_PI) phi += two_pi;
    return phi;
  }

  bool finite_point(const Point& point)
  {
    return std::isfinite(point.x) && std::isfinite(point.y) &&
           std::isfinite(point.z) && std::isfinite(point.adc);
  }

  bool solve_3x3(double A[3][3], double b[3], double x[3])
  {
    double M[3][4] = {
      {A[0][0], A[0][1], A[0][2], b[0]},
      {A[1][0], A[1][1], A[1][2], b[1]},
      {A[2][0], A[2][1], A[2][2], b[2]}
    };

    for (int col = 0; col < 3; ++col)
    {
      int pivot = col;
      for (int row = col + 1; row < 3; ++row)
      {
        if (std::fabs(M[row][col]) > std::fabs(M[pivot][col])) pivot = row;
      }
      if (std::fabs(M[pivot][col]) < 1.0e-20) return false;
      if (pivot != col)
      {
        for (int k = col; k < 4; ++k) std::swap(M[col][k], M[pivot][k]);
      }
      const double div = M[col][col];
      for (int k = col; k < 4; ++k) M[col][k] /= div;
      for (int row = 0; row < 3; ++row)
      {
        if (row == col) continue;
        const double factor = M[row][col];
        for (int k = col; k < 4; ++k) M[row][k] -= factor * M[col][k];
      }
    }

    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
  }

  double median(std::vector<double> values)
  {
    if (values.empty()) return 0.0;
    const std::size_t middle_index = values.size() / 2;
    auto middle = values.begin() + static_cast<std::ptrdiff_t>(middle_index);
    std::nth_element(values.begin(), middle, values.end());
    const double upper = *middle;
    if (values.size() % 2 != 0) return upper;
    const double lower = *std::max_element(values.begin(), middle);
    return 0.5 * (lower + upper);
  }

  double robust_sigma(const std::vector<double>& residuals, const double floor)
  {
    if (residuals.empty()) return floor;
    const double center = median(residuals);
    std::vector<double> deviations;
    deviations.reserve(residuals.size());
    for (const double residual : residuals)
    {
      deviations.push_back(std::fabs(residual - center));
    }
    return std::max(1.4826 * median(std::move(deviations)), floor);
  }

  double huber_weight(const double normalized_residual, const double threshold)
  {
    const double value = std::fabs(normalized_residual);
    if (!std::isfinite(value)) return 0.0;
    if (value <= threshold || value <= 1.0e-15) return 1.0;
    return threshold / value;
  }

  bool valid_seed(const Seed& seed)
  {
    return seed.valid && std::isfinite(seed.x) && std::isfinite(seed.y) &&
           std::isfinite(seed.z) && std::isfinite(seed.px) &&
           std::isfinite(seed.py) && std::isfinite(seed.pz) &&
           std::hypot(seed.px, seed.py) > 1.0e-12;
  }

  std::vector<std::size_t> active_indices(const std::vector<WorkPoint>& points)
  {
    std::vector<std::size_t> indices;
    indices.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      if (points[i].active) indices.push_back(i);
    }
    return indices;
  }

  double total_weight(const WorkPoint& point)
  {
    return point.base_weight * point.robust_weight;
  }

  bool algebraic_circle(const std::vector<WorkPoint>& points,
                        const std::vector<std::size_t>& indices,
                        double& xc, double& yc, double& radius)
  {
    double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double b[3] = {0.0, 0.0, 0.0};

    for (const std::size_t index : indices)
    {
      const auto& work = points[index];
      const auto& point = work.point;
      const double weight = total_weight(work);
      if (!(weight > 0.0) || !std::isfinite(weight)) continue;

      const double row[3] = {point.x, point.y, 1.0};
      const double rhs = -(point.x * point.x + point.y * point.y);
      for (int i = 0; i < 3; ++i)
      {
        b[i] += weight * row[i] * rhs;
        for (int j = 0; j < 3; ++j)
        {
          A[i][j] += weight * row[i] * row[j];
        }
      }
    }

    double solution[3] = {0.0, 0.0, 0.0};
    if (!solve_3x3(A, b, solution)) return false;

    xc = -0.5 * solution[0];
    yc = -0.5 * solution[1];
    const double radius2 = xc * xc + yc * yc - solution[2];
    if (radius2 <= 0.0) return false;
    radius = std::sqrt(radius2);
    return std::isfinite(xc) && std::isfinite(yc) &&
           std::isfinite(radius) && radius > 0.0;
  }

  void seed_circle(const Seed& seed, const double ref_xc, const double ref_yc,
                   const double ref_radius, double& xc, double& yc, double& radius)
  {
    xc = ref_xc;
    yc = ref_yc;
    radius = ref_radius;
    if (!valid_seed(seed)) return;

    const double pt = std::hypot(seed.px, seed.py);
    const double nx = -seed.py / pt;
    const double ny = seed.px / pt;
    const double cx1 = seed.x + nx * ref_radius;
    const double cy1 = seed.y + ny * ref_radius;
    const double cx2 = seed.x - nx * ref_radius;
    const double cy2 = seed.y - ny * ref_radius;
    const double d1 = std::hypot(cx1 - ref_xc, cy1 - ref_yc);
    const double d2 = std::hypot(cx2 - ref_xc, cy2 - ref_yc);
    xc = (d1 <= d2) ? cx1 : cx2;
    yc = (d1 <= d2) ? cy1 : cy2;
  }

  bool refine_circle(const std::vector<WorkPoint>& points,
                     const std::vector<std::size_t>& indices,
                     double& xc, double& yc, double& radius)
  {
    for (int iteration = 0; iteration < 20; ++iteration)
    {
      double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double b[3] = {0.0, 0.0, 0.0};
      double chi2 = 0.0;

      for (const std::size_t index : indices)
      {
        const auto& work = points[index];
        const auto& point = work.point;
        const double weight = total_weight(work);
        if (!(weight > 0.0) || !std::isfinite(weight)) continue;

        const double dx = point.x - xc;
        const double dy = point.y - yc;
        const double rho = std::hypot(dx, dy);
        if (rho <= 1.0e-12) return false;
        const double residual = rho - radius;
        const double jacobian[3] = {-dx / rho, -dy / rho, -1.0};
        for (int i = 0; i < 3; ++i)
        {
          b[i] += -weight * jacobian[i] * residual;
          for (int j = 0; j < 3; ++j)
          {
            A[i][j] += weight * jacobian[i] * jacobian[j];
          }
        }
        chi2 += weight * residual * residual;
      }

      double step[3] = {0.0, 0.0, 0.0};
      if (!solve_3x3(A, b, step)) return false;

      double scale = 1.0;
      bool accepted = false;
      for (int trial = 0; trial < 8; ++trial)
      {
        const double trial_xc = xc + scale * step[0];
        const double trial_yc = yc + scale * step[1];
        const double trial_radius = radius + scale * step[2];
        if (trial_radius > 0.0 && std::isfinite(trial_xc) &&
            std::isfinite(trial_yc) && std::isfinite(trial_radius))
        {
          double trial_chi2 = 0.0;
          for (const std::size_t index : indices)
          {
            const auto& work = points[index];
            const double weight = total_weight(work);
            const auto& point = work.point;
            const double residual =
                std::hypot(point.x - trial_xc, point.y - trial_yc) - trial_radius;
            trial_chi2 += weight * residual * residual;
          }
          if (trial_chi2 <= chi2)
          {
            xc = trial_xc;
            yc = trial_yc;
            radius = trial_radius;
            accepted = true;
            break;
          }
        }
        scale *= 0.5;
      }

      if (!accepted) return true;
      if (scale * (std::hypot(step[0], step[1]) + std::fabs(step[2])) < 1.0e-7)
      {
        return true;
      }
    }
    return true;
  }

  double direction_sign(const std::vector<WorkPoint>& points,
                        const std::vector<std::size_t>& indices,
                        const Seed& seed,
                        const double xc, const double yc, const double radius)
  {
    if (valid_seed(seed))
    {
      double best = std::numeric_limits<double>::max();
      std::size_t best_index = indices.front();
      for (const std::size_t index : indices)
      {
        const auto& point = points[index].point;
        const double distance = std::hypot(point.x - seed.x, point.y - seed.y);
        if (distance < best)
        {
          best = distance;
          best_index = index;
        }
      }
      const auto& point = points[best_index].point;
      const double rx = point.x - xc;
      const double ry = point.y - yc;
      const double tx_ccw = -ry / radius;
      const double ty_ccw = rx / radius;
      return (tx_ccw * seed.px + ty_ccw * seed.py >= 0.0) ? 1.0 : -1.0;
    }

    std::vector<double> angles;
    angles.reserve(indices.size());
    for (const std::size_t index : indices)
    {
      const auto& point = points[index].point;
      double angle = std::atan2(point.y - yc, point.x - xc);
      if (!angles.empty()) angle = unwrap_near(angle, angles.back());
      angles.push_back(angle);
    }

    double total_delta = 0.0;
    for (std::size_t i = 1; i < angles.size(); ++i)
    {
      total_delta += angles[i] - angles[i - 1];
    }
    return (total_delta >= 0.0) ? 1.0 : -1.0;
  }

  bool weighted_line_fit(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& weights,
                         double& slope, double& intercept,
                         double& chi2, int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size() || x.size() != weights.size())
    {
      return false;
    }

    double sum_w = 0.0;
    double sum_wx = 0.0;
    double sum_wy = 0.0;
    double sum_wxx = 0.0;
    double sum_wxy = 0.0;
    int used = 0;

    for (std::size_t i = 0; i < x.size(); ++i)
    {
      const double weight = weights[i];
      if (!(weight > 0.0) || !std::isfinite(weight)) continue;
      sum_w += weight;
      sum_wx += weight * x[i];
      sum_wy += weight * y[i];
      sum_wxx += weight * x[i] * x[i];
      sum_wxy += weight * x[i] * y[i];
      ++used;
    }

    const double denominator = sum_w * sum_wxx - sum_wx * sum_wx;
    if (used < 2 || sum_w <= 0.0 || std::fabs(denominator) < 1.0e-20)
    {
      return false;
    }

    slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator;
    intercept = (sum_wy - slope * sum_wx) / sum_w;
    chi2 = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
      const double residual = y[i] - (slope * x[i] + intercept);
      chi2 += weights[i] * residual * residual;
    }
    ndof = used - 2;
    return true;
  }

  bool run_fit_once(const std::vector<WorkPoint>& points,
                    const Seed& seed,
                    InternalFit& output)
  {
    output = InternalFit();
    output.active_indices = active_indices(points);
    if (output.active_indices.size() < 3) return false;

    double ref_xc = 0.0;
    double ref_yc = 0.0;
    double ref_radius = 0.0;
    if (!algebraic_circle(points, output.active_indices,
                          ref_xc, ref_yc, ref_radius)) return false;

    double xc = 0.0;
    double yc = 0.0;
    double radius = 0.0;
    seed_circle(seed, ref_xc, ref_yc, ref_radius, xc, yc, radius);
    if (!refine_circle(points, output.active_indices, xc, yc, radius)) return false;

    const double center_distance = std::hypot(xc, yc);
    if (radius <= 0.0 || center_distance <= 1.0e-12) return false;

    const double sign = direction_sign(points, output.active_indices,
                                       seed, xc, yc, radius);
    output.sign = sign;
    auto& result = output.result;
    result.xc = xc;
    result.yc = yc;
    result.radius = radius;
    result.curvature = sign / radius;
    result.d0 = sign * (center_distance - radius);

    const double perigee_x = xc * (1.0 - radius / center_distance);
    const double perigee_y = yc * (1.0 - radius / center_distance);
    const double perigee_rx = perigee_x - xc;
    const double perigee_ry = perigee_y - yc;
    result.phi0 = wrap_pi(std::atan2(sign * perigee_rx / radius,
                                    -sign * perigee_ry / radius));

    const double raw_perigee_angle = std::atan2(perigee_ry, perigee_rx);
    std::vector<double> angles;
    std::vector<double> z_values;
    std::vector<double> weights;
    angles.reserve(output.active_indices.size());
    z_values.reserve(output.active_indices.size());
    weights.reserve(output.active_indices.size());

    for (const std::size_t index : output.active_indices)
    {
      const auto& work = points[index];
      double angle = std::atan2(work.point.y - yc, work.point.x - xc);
      if (!angles.empty()) angle = unwrap_near(angle, angles.back());
      angles.push_back(angle);
      z_values.push_back(work.point.z);
      weights.push_back(total_weight(work));
    }

    double perigee_angle = unwrap_near(raw_perigee_angle, angles.front());
    while (sign * (angles.front() - perigee_angle) < 0.0)
    {
      perigee_angle -= sign * two_pi;
    }

    output.s.reserve(angles.size());
    for (const double angle : angles)
    {
      double adjusted = angle;
      while (sign * (adjusted - perigee_angle) < 0.0)
      {
        adjusted += sign * two_pi;
      }
      output.s.push_back(sign * radius * (adjusted - perigee_angle));
    }

    if (!weighted_line_fit(output.s, z_values, weights,
                           result.dzds, result.z0,
                           result.chi2_z, result.ndof_z)) return false;

    result.theta = std::atan2(1.0, result.dzds);
    result.chi2_xy = 0.0;
    output.residual_xy.reserve(output.active_indices.size());
    output.residual_z.reserve(output.active_indices.size());

    for (std::size_t i = 0; i < output.active_indices.size(); ++i)
    {
      const std::size_t index = output.active_indices[i];
      const auto& work = points[index];
      const double residual_xy =
          std::hypot(work.point.x - xc, work.point.y - yc) - radius;
      const double residual_z =
          work.point.z - (result.z0 + result.dzds * output.s[i]);
      output.residual_xy.push_back(residual_xy);
      output.residual_z.push_back(residual_z);
      result.chi2_xy += total_weight(work) * residual_xy * residual_xy;
    }

    result.ndof_xy = static_cast<int>(output.active_indices.size()) - 3;
    result.npoints_used = static_cast<int>(output.active_indices.size());
    result.ok = true;
    return true;
  }

  std::vector<WorkPoint> prepare_points(const std::vector<Point>& input,
                                        const Options& options)
  {
    std::vector<WorkPoint> points;
    points.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i)
    {
      if (finite_point(input[i])) points.push_back({input[i], i, 1.0, 1.0, true});
    }

    if (options.sort_inside_out)
    {
      std::stable_sort(points.begin(), points.end(),
          [](const WorkPoint& lhs, const WorkPoint& rhs)
          {
            return std::hypot(lhs.point.x, lhs.point.y) <
                   std::hypot(rhs.point.x, rhs.point.y);
          });
    }

    std::vector<double> positive_adc;
    positive_adc.reserve(points.size());
    for (const auto& point : points)
    {
      if (point.point.adc > 0.0) positive_adc.push_back(point.point.adc);
    }
    const double reference_adc = positive_adc.empty() ? 1.0 : median(positive_adc);

    for (auto& point : points)
    {
      if (!options.use_adc_weights)
      {
        point.base_weight = 1.0;
        continue;
      }

      if (point.point.adc <= 0.0 || reference_adc <= 0.0)
      {
        point.base_weight = options.min_adc_weight;
        continue;
      }

      const double raw_weight =
          std::pow(point.point.adc / reference_adc, options.adc_weight_power);
      point.base_weight = std::clamp(raw_weight,
                                     options.min_adc_weight,
                                     options.max_adc_weight);
    }
    return points;
  }
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points, FitResult& fit_result)
{
  Seed seed;
  Options options;
  return fit(points, seed, fit_result, options);
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points,
                             const Seed& seed,
                             FitResult& fit_result)
{
  Options options;
  return fit(points, seed, fit_result, options);
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& input,
                             const Seed& seed,
                             FitResult& fit_result,
                             const Options& options)
{
  fit_result = FitResult();
  if (options.min_points < 3) return false;

  std::vector<WorkPoint> points = prepare_points(input, options);
  if (points.size() < options.min_points) return false;

  InternalFit current;
  const int robust_iterations = std::max(0, options.max_outlier_iterations);

  for (int iteration = 0; iteration <= robust_iterations; ++iteration)
  {
    if (!run_fit_once(points, seed, current)) return false;
    if (!options.reject_outliers || iteration == robust_iterations ||
        current.active_indices.size() <= options.min_points)
    {
      break;
    }

    const double sigma_xy = robust_sigma(current.residual_xy, options.min_sigma_xy);
    const double sigma_z = robust_sigma(current.residual_z, options.min_sigma_z);

    double worst_distance = -1.0;
    std::size_t worst_local_index = 0;
    for (std::size_t i = 0; i < current.active_indices.size(); ++i)
    {
      const double normalized_xy = current.residual_xy[i] / sigma_xy;
      const double normalized_z = current.residual_z[i] / sigma_z;
      const double distance = std::hypot(normalized_xy, normalized_z);
      if (distance > worst_distance)
      {
        worst_distance = distance;
        worst_local_index = i;
      }
    }

    if (worst_distance > options.outlier_threshold &&
        current.active_indices.size() - 1 >= options.min_points)
    {
      points[current.active_indices[worst_local_index]].active = false;
      continue;
    }

    bool weights_changed = false;
    for (std::size_t i = 0; i < current.active_indices.size(); ++i)
    {
      const double normalized = std::hypot(current.residual_xy[i] / sigma_xy,
                                           current.residual_z[i] / sigma_z);
      const double new_weight = huber_weight(normalized, options.huber_threshold);
      auto& robust_weight = points[current.active_indices[i]].robust_weight;
      if (std::fabs(new_weight - robust_weight) > 1.0e-3) weights_changed = true;
      robust_weight = new_weight;
    }

    if (!weights_changed) break;
  }

  // Refit once with the final selected points and final robust weights.
  if (!run_fit_once(points, seed, current)) return false;
  fit_result = current.result;
  fit_result.npoints_rejected =
      static_cast<int>(points.size()) - fit_result.npoints_used;
  fit_result.ok = true;
  return true;
}
