#include <DGIProjection.h>
#include <DesignField.h>
#include <Logger.h>

DGIProjection::DGIProjection(Mesh &mesh_)
    : Filter(mesh_), _radius(parameters().optimization.filters.blur_filter_radius), _beta(0.1),
      _beta_max(parameters().optimization.filters.variable_thickness.dgi_projection_beta),
      _growth_rate(std::pow(
          _beta_max / _beta,
          1. / (1 + (parameters()
                         .optimization.filters.variable_thickness.dgi_projection_target_iteration -
                     parameters().optimization.continuation_start_iteration) /
                        parameters().optimization.continuation_interval))),
      _call_counter(0), _dof_handler(mesh_.tria) {}

void DGIProjection::updateContinuationParameters() {
  if (opt_iteration() < parameters().optimization.continuation_start_iteration ||
      opt_iteration() % parameters().optimization.continuation_interval != 0) {
    Logger::getInstance().addData("DGI_beta", _beta);
    return; // Only update beta at specified continuation intervals
  }
  _beta = std::min(_growth_rate * _beta, _beta_max);
  Logger::getInstance().addData("DGI_beta", _beta);
}

void DGIProjection::applyFilter(std::map<CellId, double> &pseudo_densities) {
  auto params = parameters();
  _FindNeighbors();

  // reduce input data
  std::vector<std::map<CellId, double>> temp_input =
      Utilities::MPI::all_gather(mpi_communicator(), pseudo_densities);
  std::map<CellId, double> _filtered_densities_reduced;
  for (const auto &input_map : temp_input)
    for (const auto &[cellid, val] : input_map)
      if (_filtered_densities_reduced.find(cellid) == _filtered_densities_reduced.end())
        _filtered_densities_reduced[cellid] = val;
  pseudo_densities.clear();

  pseudo_densities.clear();
  for (const auto &cell : this->_mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (std::find(this->_non_design_cells.begin(), this->_non_design_cells.end(), cell->id()) !=
        this->_non_design_cells.end()) {
      pseudo_densities[cell->id()] = _filtered_densities_reduced[cell->id()];
      continue;
    }

    _MaxDensity(cell->id(), _filtered_densities_reduced);
    _MinDensity(cell->id(), _filtered_densities_reduced);

    double diff = std::max(_max - _min, 1.e-3);
    _mean = _min + 0.5 * (_max - _min);
    _mean = std::min(std::max(_min, _mean), _max);

    _diff_map[cell->id()] = diff;
    _mean_map[cell->id()] = _mean;
    _max_map[cell->id()] = _max;
    _min_map[cell->id()] = _min;

    // Classic scaled heaviside
    double eta = (_mean - _min) / diff;
    double x = (_filtered_densities_reduced[cell->id()] - _min) / diff;
    double beta = _beta * std::pow(diff, 3);

    _beta_map[cell->id()] = beta;
    _eta_map[cell->id()] = eta;
    if (_max - _min < 1.e-3) {
      pseudo_densities[cell->id()] = _filtered_densities_reduced[cell->id()];
    } else {
      pseudo_densities[cell->id()] = diff * (tanh(beta * eta) + tanh(beta * (x - eta))) /
                                         (tanh(beta * eta) + tanh(beta * (1 - eta))) +
                                     _min;
    }

  } // end of loop over cells
  if (params.optimization.mesh_adaptivity.active &&
      _call_counter % (params.optimization.mesh_adaptivity.interval_iteration) == 0) {
    _neighbour_list.clear();
    _neighbour_list_double_r.clear();
  }
  ++_call_counter;
}

void DGIProjection::applyChainrule(std::map<CellId, double> &sensitivities) {
  _FindNeighbors();

  auto params = parameters();
  std::map<CellId, double> filtered_density_map, raw_density_map;
  filtered_density_map = DesignField::get(DesignFieldNames::FilteredDensity);

  // reduce input data
  std::vector<std::map<CellId, double>> temp_input =
      Utilities::MPI::all_gather(mpi_communicator(), sensitivities);
  std::map<CellId, double> input_reduced;
  for (const auto &input_map : temp_input)
    for (const auto &[cellid, val] : input_map)
      if (input_reduced.find(cellid) == input_reduced.end())
        input_reduced[cellid] = val;
  sensitivities.clear();

  temp_input.clear();
  temp_input = Utilities::MPI::all_gather(mpi_communicator(), filtered_density_map);
  std::map<CellId, double> filtered_densities_reduced;
  for (const auto &input_map : temp_input)
    for (const auto &[cellid, val] : input_map)
      if (filtered_densities_reduced.find(cellid) == filtered_densities_reduced.end())
        filtered_densities_reduced[cellid] = val;

  for (const auto &cell : this->_mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (std::find(this->_non_design_cells.begin(), this->_non_design_cells.end(), cell->id()) !=
        this->_non_design_cells.end()) {
      sensitivities[cell->id()] = input_reduced[cell->id()];
      continue;
    }
    _MaxDensity(cell->id(), filtered_densities_reduced);
    _MinDensity(cell->id(), filtered_densities_reduced);

    sensitivities[cell->id()] = input_reduced[cell->id()];

    if (_max - _min >= 1.e-3) {

      double diff = std::max(_max - _min, 1.e-3);
      _mean = _min + 0.5 * (_max - _min);
      _mean = std::min(std::max(_min, _mean), _max);

      // Classic scaled projection
      double eta = (_mean - _min) / diff;
      double x = (filtered_densities_reduced[cell->id()] - _min) / diff;
      double beta = _beta * std::pow(diff, 3);

      double projection_der = diff / diff * (beta - beta * std::pow(tanh(beta * (x - eta)), 2)) /
                              (tanh(beta * eta) + tanh(beta * (1 - eta)));

      sensitivities[cell->id()] *= projection_der;
    }

  } // end of loop over cells
}

void DGIProjection::_BuildCellIdCenterMap() {
  auto params = parameters();
  std::map<CellId, Point<2>> cellId_center_proc;
  std::map<CellId, double> cellId_weight_proc;
  std::map<CellId, bool> cellId_bnd_proc;
  for (const auto &cell : this->_mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    if (std::find(this->_non_design_cells.begin(), this->_non_design_cells.end(), cell->id()) !=
        this->_non_design_cells.end())
      continue;

    cellId_center_proc[cell->id()] = cell->center();
    cellId_weight_proc[cell->id()] = cell->measure();
    cellId_bnd_proc[cell->id()] = cell->at_boundary();
  }

  std::vector<std::map<CellId, Point<2>>> temp_center =
      Utilities::MPI::all_gather(mpi_communicator(), cellId_center_proc);
  std::vector<std::map<CellId, double>> temp_weight =
      Utilities::MPI::all_gather(mpi_communicator(), cellId_weight_proc);
  std::vector<std::map<CellId, bool>> temp_bnd =
      Utilities::MPI::all_gather(mpi_communicator(), cellId_bnd_proc);

  _cellId_center.clear();
  _cellId_weight.clear();
  _cellId_is_bnd.clear();

  for (const auto &center_map : temp_center)
    for (const auto &[cellid, val] : center_map)
      if (_cellId_center.find(cellid) == _cellId_center.end())
        _cellId_center[cellid] = val;
  for (const auto &weight_map : temp_weight)
    for (const auto &[cellid, val] : weight_map)
      if (_cellId_weight.find(cellid) == _cellId_weight.end())
        _cellId_weight[cellid] = val;
  for (const auto &bnd_map : temp_bnd)
    for (const auto &[cellid, val] : bnd_map)
      if (_cellId_is_bnd.find(cellid) == _cellId_is_bnd.end())
        _cellId_is_bnd[cellid] = val;
}

void DGIProjection::_FindNeighbors() {
  auto params = parameters();
  if (_neighbour_list.empty()) { // Assemble the list only once or after adaptive mesh

    _BuildCellIdCenterMap();

    for (const auto &cell : this->_mesh.tria.active_cell_iterators()) {
      if (!cell->is_locally_owned())
        continue;
      if (std::find(this->_non_design_cells.begin(), this->_non_design_cells.end(), cell->id()) !=
          this->_non_design_cells.end())
        continue;

      for (const auto &[cellid, val] : _cellId_center) {
        Point<2> center = cell->center();
        double distance = center.distance(val);
        if (distance < _radius) {
          _neighbour_list[cell->id()][cellid] = _cellId_weight[cellid];
        }
        if (distance < 2 * _radius) {
          _neighbour_list_double_r[cell->id()][cellid] = _cellId_is_bnd[cellid];
        }
      }
    }

    // Reduce neighbour list
    std::vector<std::map<CellId, std::map<CellId, double>>> temp_neighbour_list =
        Utilities::MPI::all_gather(mpi_communicator(), _neighbour_list);

    for (const auto &nlist_map : temp_neighbour_list)
      for (const auto &[cellid, val] : nlist_map)
        if (_neighbour_list.find(cellid) == _neighbour_list.end())
          _neighbour_list[cellid] = val;

    std::vector<std::map<CellId, std::map<CellId, bool>>> temp_neighbour_list_double_r =
        Utilities::MPI::all_gather(mpi_communicator(), _neighbour_list_double_r);

    for (const auto &nlist_map : temp_neighbour_list_double_r)
      for (const auto &[cellid, val] : nlist_map)
        if (_neighbour_list_double_r.find(cellid) == _neighbour_list_double_r.end())
          _neighbour_list_double_r[cellid] = val;
  }
  // end of neighbour list creation ----------------------

  // If neighbour list is still empty that means the filtering radius is smaller than the element
  // size In that case throw an error otherwise the topology optimization would produce checkboard
  // patterns
  if (_neighbour_list.empty()) {
    throw std::runtime_error(
        "Neighbour list for filtering is empty. "
        "That means the filtering radius was chosen smaller than the base element size. "
        "Try increasing the filtering radius.");
  }
}

void DGIProjection::_MaxDensity(CellId cell_id,
                                std::map<CellId, double> &reduced_filtered_densities) {
  auto params = parameters();
  _max = 0.; // params.optimization.design_field.variable_thickness.low_thickness_threshold;
  for (const auto &[n_id, value] : _neighbour_list_double_r[cell_id]) {
    if (reduced_filtered_densities.at(n_id) > _max)
      _max = reduced_filtered_densities.at(n_id);
  }
}

void DGIProjection::_MinDensity(CellId cell_id,
                                std::map<CellId, double> &reduced_filtered_densities) {
  auto params = parameters();
  _min = 1.;
  for (const auto &[n_id, value] : _neighbour_list_double_r[cell_id]) {
    if (reduced_filtered_densities.at(n_id) < _min)
      _min = reduced_filtered_densities.at(n_id);
  }
}