// C++ headers

// Deal.II headers

// Project headers
#include <DGIProjection.h>
#include <HeavisideProjection.h>
#include <Logger.h>
#include <PDEFilter.h>
#include <TopOpt.h>

using namespace dealii;

TopOpt::TopOpt()
    : _mesh(), _bvp(BVPFactory::createBVP(parameters().general.bvp_type, _mesh)),
      _optimizer(OptimizerFactory::createOptimizer()), _data_transfer(_mesh),
      _response_handler(*_bvp),
      _output(std::make_unique<DataOutput<DomainParallelTriaType, DoFHandler<2>>>("output",
                                                                                  _mesh.tria)) {}

void TopOpt::_FilterDensities() {
  TimerOutput::Scope t(compute_timer(), "TopOpt::_FilterDensities");
  auto temp_field = DesignField::get(DesignFieldNames::PseudoDensity);
  for (const auto &filter : _filters) {
    filter->applyFilter(temp_field);
    // Reset padding cells
    for (const auto &cell : _mesh.tria.active_cell_iterators()) {
      if (!cell->is_locally_owned())
        continue;
      if (cell->material_id() == e_padding_mat_id) // Padding region for hyper L geometry
        temp_field[cell->id()] = 0.0;
    }
    std::string name = filter->getFilterName();
    if (name == FilterNames::PDEFilter)
      DesignField::set(temp_field, DesignFieldNames::FilteredDensity);
    else if (name == FilterNames::DGIProjection)
      DesignField::set(temp_field, DesignFieldNames::DGIProjectedDensity);
    else if (name == FilterNames::HeavisideProjection)
      DesignField::set(temp_field, DesignFieldNames::HeavisideProjectedDensity);
  }
}

void TopOpt::_AdaptiveRefinement() {
  TimerOutput::Scope t(compute_timer(), "TopOpt::_AdaptiveRefinement");

  auto params = parameters();

  if (!(params.optimization.mesh_adaptivity.active &&
        opt_iteration() % params.optimization.mesh_adaptivity.interval_iteration == 0 &&
        opt_iteration() >= params.optimization.mesh_adaptivity.start_iteration))
    return;

  // Design fields data goes to data transfer
  auto design_fields_for_transfer = DesignField::getAllActive();
  for (auto &[name, df_data_map] : design_fields_for_transfer)
    _data_transfer.RegisterDataMap(df_data_map, name);

  // MMA data goes to data transfer
  auto mma_data_for_transfer = _optimizer->GetDataRelevantForTransfer();
  for (auto &[name, mma_data_map] : mma_data_for_transfer)
    _data_transfer.RegisterDataMap(mma_data_map, name);

  // Data transfer gets prepared
  _data_transfer.PrepareForTransfer();

  if (params.optimization.mesh_adaptivity.type == "CNF") {
    std::string cnf_str = "CNF";
    Vector<double> config_forces = _bvp->GetRefinementData(cnf_str);
    _mesh.hAdaptationCNF(config_forces, _bvp->_dof_handler);
  } else if (params.optimization.mesh_adaptivity.type == "VONMISES") {
    std::string vnm_str = "VON_MISES";
    Vector<double> von_mises = _bvp->GetRefinementData(vnm_str);
    _mesh.hAdaptationVNM(von_mises, _bvp->_smoothing_dof_handler);
  } else if (params.optimization.mesh_adaptivity.type == "DENSITY") {
    _mesh.hAdaptationDENS();
  } else if (params.optimization.mesh_adaptivity.type == "DENSITY-JUMP") {
    _mesh.hAdaptationDensityJump();
  } else {
    throw std::runtime_error("Criterion not supported");
  }

  // Data transfer execution
  _data_transfer.ExecuteTransfer();

  // Design fields receive new maps
  for (auto &[name, df_data_map] : design_fields_for_transfer)
    DesignField::set(_data_transfer.GetTransferedMap(name), name);

  // Set MMA data back to MMA class
  for (auto &[name, mma_data_map] : mma_data_for_transfer)
    _optimizer->SetDataAfterTransfer(_data_transfer.GetTransferedMap(name), name);

  // Top Opt class resets the non design cells
  _SetNonDesignCells();
}

void TopOpt::_UpdateContinuationParameters() {
  DesignField::UpdatePenalty(_change);
  Logger::getInstance().addData("Mat_penalty", DesignField::getPenalty());
  // Function applyFilter logs the filter specific continuation parameters, so we call it for all
  // filters to update the parameters and log them as well.
  for (auto &filter : _filters)
    filter->updateContinuationParameters();
}

void TopOpt::_PushMapdataToOutput(const std::map<CellId, double> &data, std::string name,
                                  DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {

  const unsigned int n_cells = _mesh.tria.n_active_cells();
  Vector<double> vec_data(n_cells);
  int i = 0;
  for (auto cell : _bvp->_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) {
      ++i;
      continue;
    }
    vec_data[i] = data.at(cell->id());
    ++i;
  }

  std::vector<std::string> names(1, name);
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(DataComponentInterpretation::component_is_scalar);
  output.template PushDataName<Vector<double>>(vec_data, names, data_component_interpretation,
                                               &_bvp->_dof_handler);
}

void TopOpt::_SetNonDesignCells() {
  auto params = parameters();
  if (_non_design_cells.empty() && !params.optimization.non_design_id.empty()) {

    std::vector<std::pair<Point<2>, Tensor<1, 2>>> boundary_face_normal_local;
    std::vector<double> non_design_thickness_local;

    FESystem<2> _fe(FE_Q<2>(params.bvp.poly_degree), 2);
    const QGauss<1> face_quadrature(_fe.degree * 2 + 1);
    FEFaceValues<2> fe_face_values(_fe, face_quadrature, update_normal_vectors);

    DoFHandler<2> dof_handler(_bvp->_tria);
    dof_handler.distribute_dofs(_fe);

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned())
        continue;

      if (cell->at_boundary()) {
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
          auto iter =
              std::find(params.optimization.non_design_id.begin(),
                        params.optimization.non_design_id.end(), cell->face(f)->boundary_id());
          if (iter != params.optimization.non_design_id.end()) {
            _non_design_cells.push_back(cell->id());

            fe_face_values.reinit(cell, f);
            boundary_face_normal_local.emplace_back(
                std::make_pair(cell->face(f)->center(), fe_face_values.normal_vector(0)));

            unsigned int idx = std::distance(params.optimization.non_design_id.begin(), iter);
            non_design_thickness_local.push_back(params.optimization.non_design_thickness[idx]);
            continue;
          }
        }
      }
    }

    // Gathering boundary_face_normal data from all the proc
    std::vector<std::pair<Point<2>, Tensor<1, 2>>> boundary_face_normal;
    std::vector<std::vector<std::pair<Point<2>, Tensor<1, 2>>>> vec_boundary_face_normal =
        Utilities::MPI::all_gather(mpi_communicator(), boundary_face_normal_local);
    for (auto &proc_boundary_face_normal : vec_boundary_face_normal)
      boundary_face_normal.insert(boundary_face_normal.end(), proc_boundary_face_normal.begin(),
                                  proc_boundary_face_normal.end());

    std::vector<double> non_design_thickness;
    std::vector<std::vector<double>> vec_non_design_thickness =
        Utilities::MPI::all_gather(mpi_communicator(), non_design_thickness_local);
    for (auto &proc_non_design_thickness : vec_non_design_thickness)
      non_design_thickness.insert(non_design_thickness.end(), proc_non_design_thickness.begin(),
                                  proc_non_design_thickness.end());

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned())
        continue;

      unsigned int idx = 0;
      for (const auto &bnd_cell : boundary_face_normal) {

        // Distance check
        if (cell->center().distance(bnd_cell.first) < non_design_thickness[idx]) {

          // Direction check
          Tensor<1, 2> distance_vector;
          distance_vector[0] = bnd_cell.first[0] - cell->center()[0];
          distance_vector[1] = bnd_cell.first[1] - cell->center()[1];
          if (std::abs(distance_vector * bnd_cell.second) >
              0.999999 * cell->center().distance(bnd_cell.first)) {
            _non_design_cells.push_back(cell->id());
          }
        } // if distance check satisified
        ++idx;
      } // end of loop over boundary cell

    } // end of loop over cells

    for (auto &filter : _filters)
      filter->SetNonDesignCells(_non_design_cells);
  }
}

void TopOpt::_SetInitialDensity() {

  auto params = parameters();

  // Initial Pseudo density value
  std::map<CellId, double> temp_dens;
  for (const auto &cell : _mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    // Stress testing setup
    // if (params.general.problem_type == "std") {
    //   // Set circular notch for L-beam test
    //   if (params.geometry.name == "hyper_L") {
    //     double notch_radius = 3.0;
    //     bool antialiasing = true;
    //     double cell_length = cell->diameter();
    //     temp_dens[cell->id()] =
    //         SetLBeamNotchDensity(cell->center(), notch_radius, antialiasing, cell_length);
    //   } else {
    //     temp_dens[cell->id()] = 1.0;
    //   }
    //   continue;
    // }

    Point<2> pt;
    pt[0] = 3;
    pt[1] = 2;

    if (std::find(_non_design_cells.begin(), _non_design_cells.end(), cell->id()) !=
        _non_design_cells.end())
      temp_dens[cell->id()] = 1.0;
    else if (cell->material_id() == e_padding_mat_id) // Padding region for hyper L geometry
      temp_dens[cell->id()] = 0.0;
    else
      temp_dens[cell->id()] = params.optimization.design_field.initial_density;
  } // end of loop over cells
  DesignField::Initialize();
  DesignField::set(temp_dens, DesignFieldNames::PseudoDensity);

} // end of _SetInitialDensity functions

void TopOpt::_SetFilters() {
  auto params = parameters();

  auto temp_dens = DesignField::get(DesignFieldNames::PseudoDensity);
  if (params.optimization.filters.blur_filter) {
    _filters.push_back(std::make_unique<PDEFilter>(_mesh, "density"));
    DesignField::set(temp_dens, DesignFieldNames::FilteredDensity);
  }

  if (params.optimization.design_field.mode == OptModeNames::VariableThickness &&
      params.optimization.filters.variable_thickness.dgi_projection) {
    _filters.push_back(std::make_unique<DGIProjection>(_mesh));
    DesignField::set(temp_dens, DesignFieldNames::DGIProjectedDensity);
  }

  if (params.optimization.filters.heaviside_projection) {
    _filters.push_back(std::make_unique<HeavisideProjection>(_mesh));
    DesignField::set(temp_dens, DesignFieldNames::HeavisideProjectedDensity);
  }
}

void TopOpt::_PrepareOptimization() {
  _SetNonDesignCells();
  _SetInitialDensity();
  _SetFilters();
  opt_iteration() = 0;
  _change = 1.0;
  PROJ_MPI_BARRIER
  Logger::getInstance().clearHistory();
  Logger::getInstance().printConsoleHeader();
  PROJ_MPI_BARRIER
}

void TopOpt::_ComputeResponseValuesGradients() {
  _response_handler.ComputeResponseValuesGradients();
}

void TopOpt::_DesignUpdate() {

  auto temp_dens = DesignField::get(DesignFieldNames::PseudoDensity);
  DesignField::set(temp_dens, DesignFieldNames::PreviousPseudoDensity);
  auto modified_values = _response_handler.GetModifiedValues();
  auto chainrule_gradients = _response_handler.GetChainruleGradients();

  _optimizer->DesignUpdate(temp_dens, chainrule_gradients, modified_values);

  for (auto &non_design_cell : this->_non_design_cells) {
    if (temp_dens.find(non_design_cell) != temp_dens.end())
      temp_dens[non_design_cell] = 1.0;
  }
  // Reset padding cells
  for (const auto &cell : _mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (cell->material_id() == e_padding_mat_id) // Padding region for hyper L geometry
      temp_dens[cell->id()] = 0.0;
  }

  DesignField::set(temp_dens, DesignFieldNames::PseudoDensity);
}

void TopOpt::_WriteDomainOutput(const unsigned int &iter) {
  PROJ_MPI_BARRIER
  _bvp->GetAllOutputData(*_output);
  if (parameters().general.problem_type == "opt") {
    _data_transfer.GetOutputData(*_output);
    _response_handler.GetAllOutputData(*_output);
  } else {
    auto all_densities = DesignField::getAllActive();
    for (auto &[name, field] : all_densities)
      _PushMapdataToOutput(field, name, *_output);
  }
  PROJ_MPI_BARRIER
  _output->WriteDataOutput(iter);
}

void TopOpt::_applyChainrule() {
  TimerOutput::Scope t(compute_timer(), "TopOpt::_applyChainrule");
  auto gradients = _response_handler.GetNormalizedGradients();
  for (auto &[func_id, grad] : gradients) {
    for (const auto &filter : _filters)
      filter->applyChainrule(grad);
  }
  _response_handler.SetChainruleGradients(gradients);
}

void TopOpt::_AddDensityDataToOutput() {

  auto params = parameters();

  auto all_densities = DesignField::getAllActive();
  for (auto &[name, field] : all_densities)
    _PushMapdataToOutput(field, name, *this->_output);

  const unsigned int n_cells = this->_mesh.tria.n_active_cells();

  auto raw_gradients = _response_handler.GetRawGradients();
  auto chainrule_gradients = _response_handler.GetChainruleGradients();

  std::map<unsigned int, Vector<double>> vec_raw_gradients, vec_chainrule_gradients;
  for (auto &[id, grad] : raw_gradients) {
    vec_raw_gradients[id].reinit(n_cells);
    vec_chainrule_gradients[id].reinit(n_cells);
  }

  unsigned int i = 0;
  for (auto &[id, grad] : raw_gradients) {
    i = 0;
    for (auto cell : _bvp->_dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) {
        ++i;
        continue;
      }
      vec_raw_gradients[id][i] = grad[cell->id()];
      ++i;
    }
  }
  for (auto &[id, grad] : chainrule_gradients) {
    i = 0;
    for (auto cell : _bvp->_dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) {
        ++i;
        continue;
      }
      vec_chainrule_gradients[id][i] = grad[cell->id()];
      ++i;
    }
  }

  std::string objective_name =
      findKeyByValue<std::string, ResponseNames>(ResponseMap, params.optimization.objective.name) +
      "_objective";
  unsigned int id = params.optimization.objective.id;

  std::vector<std::string> grad_names(1, objective_name + "_grad_raw");
  std::vector<std::string> chainrule_grad_names(1, objective_name + "_grad_chainrule");
  this->_output->template PushDataName<Vector<double>>(
      vec_raw_gradients[id], grad_names, DataComponentInterpretationTypes::scalar_interpretation,
      &_bvp->_dof_handler);
  this->_output->template PushDataName<Vector<double>>(
      vec_chainrule_gradients[id], chainrule_grad_names,
      DataComponentInterpretationTypes::scalar_interpretation, &_bvp->_dof_handler);

  for (auto constraint_params : params.optimization.constraints) {
    std::string response_name =
        findKeyByValue<std::string, ResponseNames>(ResponseMap, constraint_params.name) + "_" +
        findKeyByValue<std::string, ConstraintTypeNames>(ConstraintTypeMap, constraint_params.type);
    unsigned int id = constraint_params.id;

    std::vector<std::string> grad_names(1, response_name + "_grad_raw");
    std::vector<std::string> chainrule_grad_names(1, response_name + "_grad_chainrule");
    this->_output->template PushDataName<Vector<double>>(
        vec_raw_gradients[id], grad_names, DataComponentInterpretationTypes::scalar_interpretation,
        &_bvp->_dof_handler);
    this->_output->template PushDataName<Vector<double>>(
        vec_chainrule_gradients[id], chainrule_grad_names,
        DataComponentInterpretationTypes::scalar_interpretation, &_bvp->_dof_handler);
  }
}

void TopOpt::_WriteIterationData() {
  PROJ_MPI_BARRIER
  auto raw_values = _response_handler.GetRawValues();
  auto mod_values = _response_handler.GetModifiedValues();
  time_stamp = std::chrono::steady_clock::now(); // Assuming time_stamp is a member or local
  auto &logger = Logger::getInstance();

  logger.printConsoleDataRow(opt_iteration(),
                             raw_values, // Pass the actual map
                             mod_values, // Pass the actual map
                             _change);

  logger.addData("iteration", int(opt_iteration()));
  logger.addData("objective", raw_values[0]);
  logger.addData("change", _change);
  logger.addData("cells", int(_mesh.tria.n_global_active_cells()));
  logger.addData("dofs", int(_bvp->_dof_handler.n_dofs()));
  logger.addData(
      "time",
      double(std::chrono::duration_cast<std::chrono::seconds>(time_stamp - time_begin).count()));
  PROJ_MPI_BARRIER
}

void TopOpt::_ComputeStoppingCriterion() {
  auto params = parameters();
  double volume_per_proc = 0;
  double accum_change_per_proc = 0.0, max_change_per_proc = 0.0;
  const auto temp_dens = DesignField::get(DesignFieldNames::PseudoDensity);
  const auto temp_dens_prev = DesignField::get(DesignFieldNames::PreviousPseudoDensity);
  for (const auto &cell : _bvp->_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (temp_dens_prev.find(cell->id()) != temp_dens_prev.end()) {
      if (std::isnan(temp_dens.at(cell->id())))
        std::cerr << "Density is nan! " << std::endl;
      if (std::isnan(temp_dens_prev.at(cell->id())))
        std::cerr << "Old density is nan! " << std::endl;
      double diff = std::abs(temp_dens.at(cell->id()) - temp_dens_prev.at(cell->id()));
      accum_change_per_proc += diff * cell->measure(); // Mean criterion
      if (diff > max_change_per_proc)                  // Max criterion
        max_change_per_proc = diff;
    }
    volume_per_proc += cell->measure();
  }

  double volume = Utilities::MPI::sum(volume_per_proc, mpi_communicator());
  _change = Utilities::MPI::sum(accum_change_per_proc, mpi_communicator()) / volume;
}

void TopOpt::Run() {

  time_begin = std::chrono::steady_clock::now(); // Timer start

  auto params = parameters();

  if (params.general.problem_type == "std") {
    _SetInitialDensity(); // Set all to one since we want to run a standard analysis without
                          // optimization
    _bvp->Run();          // If we only run a BVP, then terminate afterwards
    _WriteDomainOutput(opt_iteration());
    return;
  } else if (params.general.problem_type != "opt")
    throw std::runtime_error("Wrong problem type.");

  _PrepareOptimization(); // Set non  design cells, initial densities,
                          // materials IDs, print header

  while (_change > params.optimization.convergence.design_change_tol &&
         opt_iteration() < params.optimization.convergence.max_iterations) {
    _UpdateContinuationParameters(); // Avoid updating at iteration 0 to allow logging of initial
                                     // parameters before they get updated for the first time
    _bvp->Run();
    _ComputeResponseValuesGradients();
    _applyChainrule();
    _AddDensityDataToOutput();
    _WriteDomainOutput(opt_iteration());
    _WriteIterationData();
    _DesignUpdate();
    _ComputeStoppingCriterion();
    _FilterDensities();
    _AdaptiveRefinement(); // Optional: see the function definition
    ++opt_iteration();
  }

  // Final design computation
  _UpdateContinuationParameters();
  _bvp->Run();
  _ComputeResponseValuesGradients();
  _applyChainrule();
  _AddDensityDataToOutput();
  _WriteDomainOutput(opt_iteration());
  _WriteIterationData();

  PROJ_MPI_BARRIER
  Logger::getInstance().printConsoleFooter(_change);
  Logger::getInstance().writeConvergenceHistory();
  PROJ_MPI_BARRIER
  _PrepareForPostprocessing();
}

void TopOpt::_PrepareForPostprocessing() {
  TimerOutput::Scope t(compute_timer(), "TopOpt::_PrepareForPostprocessing");

  auto params = parameters();

  for (unsigned int ref_it = 0; ref_it < params.optimization.mesh_adaptivity.max_h_refinements;
       ++ref_it) {

    // Design fields data goes to data transfer
    auto design_fields_for_transfer = DesignField::getAllActive();
    for (auto &[name, df_data_map] : design_fields_for_transfer)
      _data_transfer.RegisterDataMap(df_data_map, name);

    // MMA data goes to data transfer
    auto mma_data_for_transfer = _optimizer->GetDataRelevantForTransfer();
    for (auto &[name, mma_data_map] : mma_data_for_transfer)
      _data_transfer.RegisterDataMap(mma_data_map, name);

    // Data transfer gets prepared
    _data_transfer.PrepareForTransfer();

    // Refine all cells to eliminate hanging nodes and element size difference
    _mesh.RefineAll();

    // Data transfer execution
    _data_transfer.ExecuteTransfer();

    // Design fields receive new maps
    for (auto &[name, df_data_map] : design_fields_for_transfer)
      DesignField::set(_data_transfer.GetTransferedMap(name), name);

    // Set MMA data back to MMA class
    for (auto &[name, mma_data_map] : mma_data_for_transfer)
      _optimizer->SetDataAfterTransfer(_data_transfer.GetTransferedMap(name), name);
  }
  _AddDensityDataToOutput();
  _data_transfer.GetOutputData(*_output);
  PROJ_MPI_BARRIER
  _output->WriteDataOutput(++opt_iteration());
}