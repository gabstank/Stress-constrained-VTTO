#include <DataTransfer.h>

DataTransfer::DataTransfer(Mesh &mesh_)
    : _density_dof_handler(mesh_.tria), _fe_density(FE_DGQ<2>(0)) {}

void DataTransfer::RegisterDataMap(const std::map<CellId, double> &data_map_in,
                                   const std::string &key) {
  if (_transfer_data.find(key) == _transfer_data.end()) {
    _transfer_data.insert(
        std::make_pair(key, std::make_unique<DataUnit>(data_map_in, _density_dof_handler)));
  } else {
    _transfer_data[key]->in_out_map = data_map_in;
    _transfer_data[key]->transfered = false;
  }
}

std::map<CellId, double> DataTransfer::GetTransferedMap(const std::string &key) {
  if (_transfer_data.find(key) == _transfer_data.end())
    throw std::runtime_error("DataTransfer: For the given key there is no registered data.");
  else {
    if (!_transfer_data[key]->transfered)
      throw std::runtime_error("DataTransfer: For the given key the data has not been transfered.");
    else
      return _transfer_data[key]->in_out_map;
  }
}

void DataTransfer::PrepareForTransfer() {
  _density_dof_handler.distribute_dofs(_fe_density); // Update dofs

  _density_constraints.clear(); // Update constraints
  DoFTools::make_hanging_node_constraints(_density_dof_handler, _density_constraints);
  _density_constraints.close();

  _locally_owned_dofs = _density_dof_handler.locally_owned_dofs(); // Update locally owned dofs
  DoFTools::extract_locally_relevant_dofs(_density_dof_handler, _locally_relevant_dofs);

  for (auto &[key, data_unit] : _transfer_data) { // Reinitialize the vectors
    data_unit->in_out_vec.reinit(_locally_owned_dofs, mpi_communicator());
    data_unit->transfer_vec.reinit(_locally_owned_dofs, _locally_relevant_dofs, mpi_communicator());

    if (data_unit->transfered)
      throw std::runtime_error("DataTransfer: Following data unit have not been updated: " + key);
  }

  unsigned int dofs_per_cell = _fe_density.dofs_per_cell;
  if (dofs_per_cell != 1)
    throw std::runtime_error("DataTransfer: Density dofs are more than 1.");
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  for (const auto &cell :
       _density_dof_handler.active_cell_iterators()) { // Transfer data from map to vector
    if (!cell->is_locally_owned())
      continue;
    cell->get_dof_indices(local_dof_indices);
    for (auto &[key, data_unit] : _transfer_data)
      data_unit->in_out_vec[local_dof_indices[0]] = data_unit->in_out_map[cell->id()];
  }

  for (auto &[key, data_unit] : _transfer_data) { // Store the data in the transfer vector and
                                                  // prepare the solution transfer object
    data_unit->transfer_vec = data_unit->in_out_vec;
    data_unit->sol_transfer.prepare_for_coarsening_and_refinement(data_unit->transfer_vec);
  }
}

void DataTransfer::ExecuteTransfer() {
  _density_dof_handler.distribute_dofs(_fe_density); // Update dofs
  _locally_owned_dofs = _density_dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(_density_dof_handler, _locally_relevant_dofs);

  for (auto &[key, data_unit] : _transfer_data) { // Interpolate data
    data_unit->in_out_vec.reinit(_locally_owned_dofs, mpi_communicator());
    data_unit->in_out_map.clear();
    data_unit->sol_transfer.interpolate(data_unit->in_out_vec);
    data_unit->transfered = true;
  }

  unsigned int dofs_per_cell = _fe_density.dofs_per_cell;
  if (dofs_per_cell != 1)
    throw std::runtime_error("DataTransfer: Density dofs are more than 1.");
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  for (const auto &cell :
       _density_dof_handler.active_cell_iterators()) { // Transfer data from vector to map
    if (!cell->is_locally_owned())
      continue;
    cell->get_dof_indices(local_dof_indices);
    for (auto &[key, data_unit] : _transfer_data)
      data_unit->in_out_map[cell->id()] = data_unit->in_out_vec[local_dof_indices[0]];
  }
}

void DataTransfer::GetOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {

  std::vector<std::string> density_name(1, "nodal_dens");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      scalar_data_component_interpretation(DataComponentInterpretation::component_is_scalar);
  if (_transfer_data.find(DesignFieldNames::PseudoDensity) != _transfer_data.end()) {
    output.template PushDataName<TrilinosWrappers::MPI::Vector>(
        _transfer_data[DesignFieldNames::PseudoDensity]->in_out_vec, density_name,
        scalar_data_component_interpretation, &_density_dof_handler);
  }

} // end of function