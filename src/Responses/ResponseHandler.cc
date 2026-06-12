// C++ headers

// Deal.II headers

// Project headers
#include <ResponseHandler.h>
#include <Volume.h>

// Linear elastic
#include <ComplianceLinear.h>
#include <PNStressLinear.h>

using namespace dealii;

ResponseHandler::ResponseHandler(BVP &bvp_) : _bvp(bvp_) {
  _AddResponse(parameters().optimization.objective.id,
               parameters().optimization.objective.name); // objective
  for (auto constr_params : parameters().optimization.constraints)
    _AddResponse(constr_params.id, constr_params.name);
}

void ResponseHandler::_ComputeValues() {

  TimerOutput::Scope t(compute_timer(), "ResponseHandler::ComputeValues");

  auto params = parameters();

  _raw_values.clear();

  for (auto &[id, response] : _responses)
    _raw_values[id] = response->GetFunction();
  if (!_reference_computed) {
    _ComputeRefValues();
    _reference_computed = true;
  }
}

void ResponseHandler::_ComputeRefValues() {
  auto params = parameters();

  _ref_values.clear();

  for (auto &[id, response] : _responses)
    _ref_values[id] = response->GetFunctionRef();
}

void ResponseHandler::_NormalizeValuesGradients() {
  auto params = parameters();

  _mod_values.clear();
  _mod_values = _raw_values;
  _normalized_gradients.clear();
  _normalized_gradients = _raw_gradients;

  // normalizing wrt to reference (initial) response values considering type of constraint
  for (auto &[id, mod_value] : _mod_values)
    mod_value /= std::abs(_ref_values[id]);
  for (auto &[id, sens] : _normalized_gradients)
    for (auto &[cellid, grad_val] : sens)
      grad_val /= std::abs(_ref_values[id]);

  // Adjusting for types of constraints
  for (auto constr_params : params.optimization.constraints) {
    if (constr_params.type == ConstraintTypeNames::equality ||
        constr_params.type == ConstraintTypeNames::upper)
      _mod_values[constr_params.id] -= 1.;
    else if (constr_params.type == ConstraintTypeNames::lower)
      _mod_values[constr_params.id] = 1. - _mod_values[constr_params.id];
  }
  for (auto constr_params : params.optimization.constraints)
    if (constr_params.type == ConstraintTypeNames::lower)
      for (auto &[cellid, grad_val] : _normalized_gradients[constr_params.id])
        grad_val *= -1;

  if (params.optimization.optimizer == "MMA")
    for (auto &[func_id, sens] : _normalized_gradients)
      L2NormalizeMap(sens);
}

void ResponseHandler::_ComputeGradients() {

  TimerOutput::Scope t(compute_timer(), "ResponseHandler::ComputeGradients");

  auto params = parameters();

  // Store the sensitivities in the map _density_gradients: response id -> map
  // of gradient values
  _raw_gradients.clear();
  double cell_gradient_density;

  // Iterate over cells
  for (const auto &cell : _bvp._dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    /*if (cell->material_id() == e_padding_mat_id) { // Skip cells in the padding region
      _raw_gradients[params.optimization.objective.id][cell->id()] = 0.;
      for (auto constr_params : params.optimization.constraints)
        _raw_gradients[constr_params.id][cell->id()] = 0.;
      continue;
    }*/ // might cause problems
    _responses[params.optimization.objective.id]->ComputeCellDensityGradient(cell,
                                                                             cell_gradient_density);
    _raw_gradients[params.optimization.objective.id][cell->id()] =
        cell_gradient_density / cell->measure();
    for (auto constr_params : params.optimization.constraints) {
      _responses[constr_params.id]->ComputeCellDensityGradient(cell, cell_gradient_density);
      _raw_gradients[constr_params.id][cell->id()] = cell_gradient_density / cell->measure();
    }
  }
}

void ResponseHandler::_AddResponse(const unsigned int &id, const ResponseNames &response_name) {
  auto params = parameters();
  if (response_name == ResponseNames::compliance) {
    if (params.general.bvp_type == "elasticityLin")
      _responses[id] = std::make_unique<ComplianceLinear>(_bvp, id);
  } else if (response_name == ResponseNames::volume) {
    _responses[id] = std::make_unique<Volume>(_bvp, id);
  } else if (response_name == ResponseNames::PNStress) {
    if (params.general.bvp_type == "elasticityLin")
      _responses[id] = std::make_unique<PNStressLinear>(_bvp, id);
  } else
    throw std::runtime_error("Runtime error @ ResponseHandler::_AddResponse:- Wrong response "
                             "name.");
  _responses_names[id] = findKeyByValue<std::string, ResponseNames>(ResponseMap, response_name);
  _response_id.push_back(id);
}

void ResponseHandler::ComputeResponseValuesGradients() {
  _ComputeValues();
  _ComputeGradients();
  _NormalizeValuesGradients();
}