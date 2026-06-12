#include <DesignField.h>
#include <Logger.h>
#include <MacrosAndTypedefs.h>
#include <Parameter.h>
#include <VariableThickness.h>

// Definition and initialization of static member variables
std::unique_ptr<Penalization> DesignField::_model = nullptr;
double DesignField::_penalty = 0;
double DesignField::_stress_penalty = 0;
std::map<CellId, double> DesignField::_pseudo_densities;
std::map<CellId, double> DesignField::_filtered_densities;
std::map<CellId, double> DesignField::_vt_projected_densities;
std::map<CellId, double> DesignField::_heaviside_projected_densities;
std::map<CellId, double> DesignField::_prev_pseudo_densities;

void DesignField::Initialize() {
  auto mode = parameters().optimization.design_field.mode;

  if (mode == OptModeNames::VariableThickness) {
    _model = std::make_unique<VariableThickness>();
  } else {
    throw std::runtime_error("Unknown Optimization Mode");
  }
  _model->InitializePenalty(_penalty);
  _model->InitializeStressExponent(_stress_penalty);
}

void DesignField::set(const std::map<CellId, double> &design_field, const std::string &name) {
  if (name == DesignFieldNames::PseudoDensity) {
    _pseudo_densities.clear();
    _pseudo_densities = design_field;
  } else if (name == DesignFieldNames::FilteredDensity) {
    _filtered_densities.clear();
    _filtered_densities = design_field;
  } else if (name == DesignFieldNames::DGIProjectedDensity) {
    _vt_projected_densities.clear();
    _vt_projected_densities = design_field;
  } else if (name == DesignFieldNames::HeavisideProjectedDensity) {
    _heaviside_projected_densities.clear();
    _heaviside_projected_densities = design_field;
  } else if (name == DesignFieldNames::PreviousPseudoDensity) {
    _prev_pseudo_densities.clear();
    _prev_pseudo_densities = design_field;
  } else
    throw std::runtime_error("Wrong design field name given to DesignField::set");
}

std::map<CellId, double> DesignField::get(const std::string &name) {
  if (name == DesignFieldNames::PseudoDensity) {
    if (_pseudo_densities.empty())
      throw std::runtime_error("PseudoDensity is empty!");
    return _pseudo_densities;
  } else if (name == DesignFieldNames::FilteredDensity) {
    if (_filtered_densities.empty())
      throw std::runtime_error("FilteredDensity is empty!");
    return _filtered_densities;
  } else if (name == DesignFieldNames::DGIProjectedDensity) {
    if (_vt_projected_densities.empty())
      throw std::runtime_error("DGIProjectedDensity is empty!");
    return _vt_projected_densities;
  } else if (name == DesignFieldNames::HeavisideProjectedDensity) {
    if (_heaviside_projected_densities.empty())
      throw std::runtime_error("HeavisideProjectedDensity is empty!");
    return _heaviside_projected_densities;
  } else if (name == DesignFieldNames::PreviousPseudoDensity) {
    if (_prev_pseudo_densities.empty())
      throw std::runtime_error("PreviousPseudoDensity is empty!");
    return _prev_pseudo_densities;
  } else
    throw std::runtime_error("Wrong design field name requested from DesignField::get");
}

std::map<CellId, double> DesignField::getCurrent() {
  if (!_heaviside_projected_densities.empty())
    return _heaviside_projected_densities;
  else if (!_vt_projected_densities.empty())
    return _vt_projected_densities;
  else if (!_filtered_densities.empty())
    return _filtered_densities;
  else if (!_pseudo_densities.empty())
    return _pseudo_densities;
  else
    throw std::runtime_error("All density maps are empty!");
}

std::map<std::string, std::map<CellId, double>> DesignField::getAllActive() {
  std::map<std::string, std::map<CellId, double>> all_active_fields;
  if (!_pseudo_densities.empty())
    all_active_fields[DesignFieldNames::PseudoDensity] = _pseudo_densities;
  if (!_filtered_densities.empty())
    all_active_fields[DesignFieldNames::FilteredDensity] = _filtered_densities;
  if (!_vt_projected_densities.empty())
    all_active_fields[DesignFieldNames::DGIProjectedDensity] = _vt_projected_densities;
  if (!_heaviside_projected_densities.empty())
    all_active_fields[DesignFieldNames::HeavisideProjectedDensity] = _heaviside_projected_densities;
  if (!_prev_pseudo_densities.empty())
    all_active_fields[DesignFieldNames::PreviousPseudoDensity] = _prev_pseudo_densities;
  return all_active_fields;
}

void DesignField::UpdatePenalty(const double design_change) {
  auto params = parameters();
  if (!_model)
    throw std::runtime_error("DesignField not initialized!");
  // Pass global parameters if needed, or if strategy stores them
  _model->UpdatePenalty(_penalty, design_change);
  _model->UpdateStressExponent(_stress_penalty);
}

double DesignField::GetPhysicalDensity(const double density) {
  return _model->GetPhysicalDensity(density, _penalty);
}

double DesignField::GetPhysicalDensityDer(const double density) {
  return _model->GetPhysicalDensityDer(density, _penalty);
}

double DesignField::StressRelaxation(const double density) {
  return _model->GetStressRelaxation(density, _stress_penalty);
}

double DesignField::StressRelaxationDer(const double density) {
  return _model->GetStressRelaxationDer(density, _stress_penalty);
}

bool DesignField::isSolidCell(const double density) { return _model->isSolidCell(density); }
