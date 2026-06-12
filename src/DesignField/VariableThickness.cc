#include <Logger.h>
#include <Parameter.h>
#include <VariableThickness.h>

void VariableThickness::InitializePenalty(double &current_penalty) {
  auto params = parameters();
  if (!params.optimization.design_field.variable_thickness.penalty_continuation) {
    current_penalty = params.optimization.design_field.penalty;
    return;
  }
  if (current_penalty < 1.e-10)
    current_penalty = 1.0;
  unsigned int n_update_steps =
      1 + (params.optimization.design_field.variable_thickness.penalty_target_iteration -
           params.optimization.continuation_start_iteration) /
              params.optimization.continuation_interval;
  _penalty_growth_rate =
      std::pow(params.optimization.design_field.penalty / current_penalty, 1.0 / n_update_steps);
}

void VariableThickness::UpdatePenalty(double &current_penalty, const double) {
  auto params = parameters();
  if (!params.optimization.design_field.variable_thickness.penalty_continuation) {
    return; // If no penalty continuation, set penalty to final value and return
  }
  if (opt_iteration() < params.optimization.continuation_start_iteration ||
      opt_iteration() % params.optimization.continuation_interval != 0)
    return; // Only update penalty at specified continuation intervals
  current_penalty *= _penalty_growth_rate;
  if (current_penalty > params.optimization.design_field.penalty)
    current_penalty = params.optimization.design_field.penalty;
}

void VariableThickness::InitializeStressExponent(double &current_exponent) {
  auto params = parameters();
  if (current_exponent < 1.e-10)
    current_exponent = 1.0;
  unsigned int n_update_steps =
      1 + (params.optimization.design_field.stress_relaxation_target_iteration -
           params.optimization.continuation_start_iteration) /
              params.optimization.continuation_interval;
  _target_stress_relaxation_exponent =
      std::log(params.optimization.design_field.stress_relaxation_epsilon_target) /
      std::log(params.optimization.design_field.stress_relaxation_epsilon_initial);
  _stress_relaxation_exponent_growth_rate =
      std::pow(_target_stress_relaxation_exponent / current_exponent, 1.0 / n_update_steps);
}

void VariableThickness::UpdateStressExponent(double &current_exponent) {
  auto params = parameters();
  if (opt_iteration() < params.optimization.continuation_start_iteration ||
      opt_iteration() % params.optimization.continuation_interval != 0) {
    Logger::getInstance().addData("Stress_exponent", current_exponent);
    Logger::getInstance().addData(
        "Stress_epsilon",
        std::pow(params.optimization.design_field.stress_relaxation_epsilon_initial,
                 current_exponent));
    return; // Only update stress relaxation exponent at specified continuation intervals
  }
  current_exponent *= _stress_relaxation_exponent_growth_rate;
  if (current_exponent > _target_stress_relaxation_exponent)
    current_exponent = _target_stress_relaxation_exponent;
  Logger::getInstance().addData("Stress_exponent", current_exponent);
  Logger::getInstance().addData(
      "Stress_epsilon", std::pow(params.optimization.design_field.stress_relaxation_epsilon_initial,
                                 current_exponent));
}

double VariableThickness::GetPhysicalDensity(const double density, const double penalty) {
  double low_thickness =
      parameters().optimization.design_field.variable_thickness.low_thickness_threshold;
  if (density >= low_thickness)
    return density;
  else {
    double local_density = density / low_thickness;
    double local_simp = std::pow(local_density, penalty);
    return local_simp * low_thickness * (1. - 1e-9) + 1e-9;
  }
}

double VariableThickness::GetPhysicalDensityDer(const double density, const double penalty) {
  double low_thickness =
      parameters().optimization.design_field.variable_thickness.low_thickness_threshold;
  if (density >= low_thickness)
    return 1;
  else {
    double local_density = density / low_thickness;
    double local_simp_der = penalty * std::pow(local_density, penalty - 1.);
    return local_simp_der * (1. - 1e-9);
  }
}

double VariableThickness::GetStressRelaxation(const double density, const double penalty) {
  auto params = parameters();
  double epsilon = params.optimization.design_field.stress_relaxation_epsilon_initial;
  epsilon = std::pow(epsilon, penalty);

  // New strategy: hybrid function that has a epsilon relaxation for solid densities and
  // linear function for low densities
  double low_thickness =
      params.optimization.design_field.variable_thickness.low_thickness_threshold;
  double smoothness_param =
      epsilon * epsilon; // This parameter controls how smooth the transition is between the linear
                         // and shifted epsilon relaxation, and should be small enough to not affect
                         // the results but large enough to avoid numerical issues with the square
                         // root. Setting it equal to epsilon^2 seems to work well in practice.
  double smooth_maximum =
      0.5 * (density + low_thickness +
             std::sqrt(std::pow(density - low_thickness, 2) + smoothness_param));

  double stress_relaxation = density / (epsilon * (1. - density) + smooth_maximum);
  return stress_relaxation;
}

double VariableThickness::GetStressRelaxationDer(const double density, const double penalty) {
  auto params = parameters();
  double epsilon = params.optimization.design_field.stress_relaxation_epsilon_initial;
  epsilon = std::pow(epsilon, penalty);

  // Derivative of the hybrid function, computed using the quotient rule and chain rule
  double low_thickness =
      params.optimization.design_field.variable_thickness.low_thickness_threshold;
  double smoothness_param = epsilon * epsilon;
  double sqrt_term = std::sqrt(std::pow(density - low_thickness, 2) + smoothness_param);
  double smooth_maximum = 0.5 * (density + low_thickness + sqrt_term);
  double dsmooth_maximum_ddensity = 0.5 * (1 + (density - low_thickness) / sqrt_term);

  double stress_relaxation_der = (epsilon * (1. - density) + smooth_maximum -
                                  density * (-epsilon + dsmooth_maximum_ddensity)) /
                                 std::pow(epsilon * (1. - density) + smooth_maximum, 2);
  return stress_relaxation_der;
}

bool VariableThickness::isSolidCell(const double density) {
  auto params = parameters();
  double low_thickness =
      params.optimization.design_field.variable_thickness.low_thickness_threshold;
  return density >= low_thickness;
}