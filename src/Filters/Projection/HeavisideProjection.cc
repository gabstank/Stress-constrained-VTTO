#include <DesignField.h>
#include <HeavisideProjection.h>
#include <Logger.h>
#include <VarThicknessHeaviside.h>

HeavisideProjection::HeavisideProjection(Mesh &mesh)
    : Filter(mesh), _beta(1),
      _beta_max(parameters().optimization.filters.heaviside_projection_beta),
      _growth_rate(std::pow(
          _beta_max / _beta,
          1. / (1 + (parameters().optimization.filters.heaviside_projection_target_iteration -
                     parameters().optimization.continuation_start_iteration) /
                        parameters().optimization.continuation_interval))) {
  auto mode = parameters().optimization.design_field.mode;

  if (mode == OptModeNames::VariableThickness) {
    _kernel = std::make_unique<VarThicknessHeaviside>();
  } else
    throw std::runtime_error("Wrong optimization mode.");
}

void HeavisideProjection::updateContinuationParameters() {
  if (opt_iteration() < parameters().optimization.continuation_start_iteration ||
      opt_iteration() % parameters().optimization.continuation_interval != 0) {
    Logger::getInstance().addData("H_beta", _beta);
    return; // Only update beta at specified continuation intervals
  }
  _beta = std::min(_growth_rate * _beta, _beta_max);
  Logger::getInstance().addData("H_beta", _beta);
}

void HeavisideProjection::applyFilter(std::map<CellId, double> &density_map) {
  _in_out_map.clear();
  _in_out_map = density_map;
  density_map.clear();

  // Delegate math to the kernel
  _kernel->Apply(density_map, _in_out_map, _beta, this->_non_design_cells, this->_mesh);
}

void HeavisideProjection::applyChainrule(std::map<CellId, double> &sensitivity_map) {
  _in_out_map.clear();
  _in_out_map = sensitivity_map;
  sensitivity_map.clear();
  std::map<CellId, double> density_map = DesignField::get(DesignFieldNames::FilteredDensity);

  // Delegate math to the kernel
  _kernel->ChainRule(sensitivity_map, _in_out_map, density_map, _beta, this->_non_design_cells,
                     this->_mesh);
}