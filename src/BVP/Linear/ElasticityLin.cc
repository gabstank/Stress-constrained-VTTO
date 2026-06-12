// C++ headers

// Deal.II headers

// Project headers
#include <DesignField.h>
#include <ElasticityLin.h>

ElasticityLin::ElasticityLin(Mesh &mesh_) : BVP(mesh_), _mat() {}

ElasticityLin::~ElasticityLin() {
  this->_dof_handler.clear();
  this->_smoothing_dof_handler.clear();
}

void ElasticityLin::_SetupSystem() {
  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_SetupSystem");

  auto params = parameters();

  this->_dof_handler.distribute_dofs(this->_fe_collection);

  this->_locally_owned_dofs = this->_dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->_dof_handler, this->_locally_relevant_dofs);

  this->_solution.reinit(this->_dof_handler.n_dofs());
  this->_system_rhs.reinit(this->_locally_owned_dofs, mpi_communicator());
  this->_external_force.reinit(this->_locally_owned_dofs, mpi_communicator());

  _SetupConstraints();
  DynamicSparsityPattern dsp(this->_locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->_dof_handler, dsp, this->_constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp, this->_dof_handler.locally_owned_dofs(),
                                             mpi_communicator(), this->_locally_relevant_dofs);
  this->_system_matrix.reinit(this->_locally_owned_dofs, this->_locally_owned_dofs, dsp,
                              mpi_communicator());

  this->_SetupSmoothingSystem();
  _SetupQPointHistory();

} // end of _SetupSystem() function

void ElasticityLin::_SetupConstraints() {
  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_SetupConstraints");

  auto params = parameters();

  this->_constraints.clear();
  this->_constraints.reinit(this->_locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->_dof_handler, this->_constraints);

  for (unsigned int i = 0; i < params.bvp.dirichlet.id.size(); ++i) {
    std::vector<bool> comp_vec(2);
    comp_vec[params.bvp.dirichlet.component[i]] = true;
    ComponentMask comp_mask(comp_vec);
    VectorTools::interpolate_boundary_values(
        this->_dof_handler, params.bvp.dirichlet.id[i],
        Functions::ConstantFunction<2>(params.bvp.dirichlet.value[i], 2), this->_constraints,
        comp_mask);
  } // loop over all the given dbc id

  this->_constraints.close();
} // end of _SetupConstraints()

void ElasticityLin::_SetupQPointHistory() {
  std::vector<QPointHistory> tmp;
  tmp.swap(quadrature_point_history);

  quadrature_point_history.resize(this->_tria.n_active_cells() *
                                  this->_q_collection.max_n_quadrature_points());

  hp::FEValues<2> fe_values_hp(this->_fe_collection, this->_q_collection,
                               update_values | update_quadrature_points);

  unsigned int history_index = 0;
  for (auto &cell : this->_tria.active_cell_iterators()) {
    fe_values_hp.reinit(cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    if (!cell->is_locally_owned()) {
      history_index += fe_values.n_quadrature_points;
      continue;
    }
    cell->set_user_pointer(&quadrature_point_history[history_index]);
    history_index += fe_values.n_quadrature_points;
  }
}

void ElasticityLin::_UpdateQPointHistory() {

  auto params = parameters();

  hp::FEValues<2> fe_values_hp(this->_fe_collection, this->_q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

  const FEValuesExtractors::Vector u_fe(0);

  auto densities = DesignField::getCurrent();

  for (auto &cell : this->_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    double cell_density = 1.;
    if (params.general.problem_type == "opt")
      cell_density = densities[cell->id()];

    auto *local_quadrature_points_history = reinterpret_cast<QPointHistory *>(cell->user_pointer());

    fe_values_hp.reinit(cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    std::vector<SymmetricTensor<2, 2>> local_displacement_gradients(fe_values.n_quadrature_points);
    fe_values[u_fe].get_function_symmetric_gradients(this->_solution, local_displacement_gradients);

    double cell_relaxation = 1.;
    if (params.general.problem_type == "opt")
      cell_relaxation = DesignField::StressRelaxation(cell_density);

    for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
      SymmetricTensor<2, 2> sym_grad_U_ = local_displacement_gradients[q];
      local_quadrature_points_history[q].sym_grad_U = sym_grad_U_;
      local_quadrature_points_history[q].cauchy = _mat.GetCauchyStress(sym_grad_U_, cell_density);
      local_quadrature_points_history[q].vm_stress =
          _mat.GetVonMisesStress(sym_grad_U_, cell_density);
      local_quadrature_points_history[q].cauchy =
          cell_relaxation * _mat.GetCauchyStress(sym_grad_U_);
      local_quadrature_points_history[q].vm_stress =
          cell_relaxation * _mat.GetVonMisesStress(sym_grad_U_);
    }
  }
}

void ElasticityLin::_AssembleSystem() {
  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_AssembleSystem()");

  auto params = parameters();

  this->_system_matrix = 0.0;
  this->_system_rhs = 0.0;

  hp::FEValues<2> fe_values_hp(this->_fe_collection, this->_q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

  const FEValuesExtractors::Vector u_fe(0);

  unsigned int dofs_per_cell, n_q_points;

  FullMatrix<double> cell_matrix;
  std::vector<types::global_dof_index> local_dof_indices;

  SymmetricTensor<4, 2> elasticity_tensor;

  auto densities = DesignField::getCurrent();

  for (const auto &cell : this->_dof_handler.active_cell_iterators()) {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;

    elasticity_tensor = _mat.GetElasticityTensor();
    double physical_density = DesignField::GetPhysicalDensity(densities[cell->id()]);

    cell_matrix = 0.0;

    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    n_q_points = fe_values.n_quadrature_points;
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        SymmetricTensor<2, 2> Grad_N_i = symmetrize(fe_values[u_fe].gradient(i, q));
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          SymmetricTensor<2, 2> Grad_N_j = symmetrize(fe_values[u_fe].gradient(j, q));
          cell_matrix(i, j) += contract3(Grad_N_i, elasticity_tensor, Grad_N_j) * fe_values.JxW(q) *
                               physical_density;
        } // end of j loop
      } // end of i loop
    } // end of loop over quadrature points

    cell->get_dof_indices(local_dof_indices);

    Vector<double> dummy_cell_rhs(dofs_per_cell);

    this->_constraints.distribute_local_to_global(cell_matrix, dummy_cell_rhs, local_dof_indices,
                                                  this->_system_matrix, this->_system_rhs);

  } // end of loop over cells

  this->_system_matrix.compress(VectorOperation::add);
  this->_system_rhs.compress(VectorOperation::add);

} // end of _AssembleSystem function

void ElasticityLin::_AssembleRhs() {
  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_AssembleRhs");

  auto params = parameters();

  this->_external_force = 0.0;

  // if no neumann boundary condition, do nothing.
  if (params.bvp.neumann.id.size() == 0)
    return;

  hp::FEFaceValues<2> fe_face_values_hp(this->_fe_collection, this->_face_q_collection,
                                        update_values | update_quadrature_points |
                                            update_JxW_values);
  hp::FEValues<2> fe_values_hp(this->_fe_collection, this->_q_collection,
                               update_values | update_quadrature_points);

  unsigned int dofs_per_cell, n_q_points;

  std::vector<unsigned int> local_dof_indices;
  Vector<double> cell_rhs;

  unsigned int neumann_id = 0;
  int neumann_comp = 0;

  // create a list of neumann_id with no repetition of neumann id.
  std::vector<int> unique_neumann_id = params.bvp.neumann.id;
  sort(unique_neumann_id.begin(), unique_neumann_id.end());
  unique_neumann_id.erase(std::unique(unique_neumann_id.begin(), unique_neumann_id.end()),
                          unique_neumann_id.end());

  // create a map with neumann id vs load area
  std::map<int, double> neumann_id_area;
  for (auto id : unique_neumann_id)
    neumann_id_area[id] = 0.0;

  for (const auto &cell : this->_dof_handler.active_cell_iterators()) {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;
    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);

    if (cell->at_boundary()) {
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        fe_face_values_hp.reinit(cell, f);
        const FEFaceValues<2> &fe_face_values = fe_face_values_hp.get_present_fe_values();
        n_q_points = fe_face_values.n_quadrature_points;
        for (auto id : unique_neumann_id) {
          if (cell->face(f)->boundary_id() == (unsigned int)id) {
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
              neumann_id_area[id] += fe_face_values.JxW(q_point);
            }
          }
        }
      }
    }
  }

  // create a map with id vs load area
  std::map<int, double> total_neumann_id_area;
  for (auto id : unique_neumann_id) {
    total_neumann_id_area[id] = Utilities::MPI::sum(neumann_id_area[id], mpi_communicator());
  }

  // actual run to assemble rhs vector
  for (const auto &cell : this->_dof_handler.active_cell_iterators()) {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;
    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    cell_rhs.reinit(dofs_per_cell);

    if (cell->at_boundary()) {
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        fe_face_values_hp.reinit(cell, f);
        const FEFaceValues<2> &fe_face_values = fe_face_values_hp.get_present_fe_values();
        n_q_points = fe_face_values.n_quadrature_points;
        for (unsigned int id = 0; id < params.bvp.neumann.id.size(); ++id) {
          neumann_id = params.bvp.neumann.id[id];
          neumann_comp = params.bvp.neumann.component[id];
          if (cell->face(f)->boundary_id() == neumann_id) {
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
              for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                const int component_i = cell->get_fe().system_to_component_index(i).first;
                if (component_i == neumann_comp)
                  cell_rhs(i) +=
                      params.bvp.neumann.value[id] * fe_face_values.shape_value(i, q_point) *
                      fe_face_values.JxW(q_point); // / total_neumann_id_area[neumann_id];
              }
            }
          }
        }
      }
    }

    cell->get_dof_indices(local_dof_indices);

    // we need to change this in order to allow for inhomogeneous dirichlet bc
    // constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
    // system_rhs);
    this->_constraints.distribute_local_to_global(cell_rhs, local_dof_indices, this->_system_rhs);
    this->_constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                                  this->_external_force);
  }

  this->_system_rhs.compress(VectorOperation::add);
  this->_external_force.compress(VectorOperation::add);
}

void ElasticityLin::_SolveSystem() {

  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_SolveSystem");

  this->_SolveLinearSystem(this->_system_matrix, this->_solution, this->_system_rhs);
}

void ElasticityLin::_Postprocess() {

  TimerOutput::Scope t(compute_timer(), "ElasticityLin::_Postprocess");

  // Standard postprocessing
  _UpdateQPointHistory();
  this->AssembleSmoothingMatrix();
  unsigned int i = e_von_mises_s;
  _AssembleSmoothingRhs(i);
  this->SolveSmoothingSystem();
  _UpdatePostprocessingData(i);

  // Elemental stress postprocessing
  _SetElementStresses();
  this->_el_stress_map_smoothed.clear();
  this->_el_stress_map_smoothed = this->_el_stress_map_raw;
  _stress_filter.applyFilter(this->_el_stress_map_smoothed);
}

void ElasticityLin::GetAllOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {
  BVP::GetAllOutputData(output);

  auto di_scalar = DataComponentInterpretationTypes::scalar_interpretation;
  std::vector<std::string> vm_names(1, "nodal_vm_stress");
  output.template PushDataName<Vector<double>>(postprocessing_data.domain_vm_stress, vm_names,
                                               di_scalar, &this->_smoothing_dof_handler);
}

Vector<double> ElasticityLin::GetRefinementData(const std::string &ref_type) {
  if (ref_type == "VON_MISES") {
    return postprocessing_data.domain_vm_stress;
  } else {
    Vector<double> vec;
    throw std::runtime_error("Wrong refinement data request.");
    return vec;
  }
}

void ElasticityLin::_SetElementStresses() {
  _el_stress_map_raw.clear();
  auto densities = DesignField::getCurrent();
  hp::QCollection<2> q_collection_midpoint;
  q_collection_midpoint.push_back(QMidpoint<2>());
  hp::FEValues<2> fe_values_hp(_fe_collection, q_collection_midpoint,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
  const FEValuesExtractors::Vector u_fe(0);
  for (auto &cell : _dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    double cell_density = densities.at(cell->id());
    double cell_stress_relaxation = DesignField::StressRelaxation(cell_density);

    fe_values_hp.reinit(cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    unsigned int n_q_points = fe_values.n_quadrature_points;

    std::vector<SymmetricTensor<2, 2>> solution_sym_grads_u(n_q_points);
    fe_values[u_fe].get_function_symmetric_gradients(this->_solution, solution_sym_grads_u);

    for (unsigned int q = 0; q < n_q_points; ++q) {
      double current_vm_stress = _mat.GetVonMisesStress(solution_sym_grads_u[q]);

      double relaxed_stress = cell_stress_relaxation * current_vm_stress;
      //  We have only one q point so we can increment the cell counter here for the output
      _el_stress_map_raw[cell->id()] = relaxed_stress;
    } // end of loop over quadrature points
  } // loop over cells
}

void ElasticityLin::_AssembleSmoothingRhs(unsigned int assembly_flag) {
  this->_smoothing_rhs = 0;

  hp::FEValues<2> fe_values_hp(this->_smoothing_fe_collection, this->_smoothing_q_collection,
                               update_values | update_quadrature_points | update_JxW_values);

  Vector<double> cell_rhs;
  std::vector<unsigned int> local_dof_indices;
  double q_value;

  for (auto &cell : this->_smoothing_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    fe_values_hp.reinit(cell);
    cell_rhs.reinit(cell->get_fe().dofs_per_cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();

    auto *local_quadrature_points_history = reinterpret_cast<QPointHistory *>(cell->user_pointer());

    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

    for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
      q_value = 0;

      if (assembly_flag == e_von_mises_s)
        q_value = local_quadrature_points_history[q].vm_stress;

      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        cell_rhs(i) += q_value * fe_values.shape_value(i, q) * fe_values.JxW(q);
      }
    }
    // add to system smoothing rhs
    local_dof_indices.resize(cell->get_fe().dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
    this->_smoothing_constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                                            this->_smoothing_rhs);
  }
  this->_smoothing_rhs.compress(VectorOperation::add);
}

void ElasticityLin::_UpdatePostprocessingData(unsigned int assembly_flag) {

  Vector<double> reduced_smoothing_solution(this->_smoothing_solution);
  if (assembly_flag == e_von_mises_s) {
    postprocessing_data.domain_vm_stress.reinit(this->_smoothing_dof_handler.n_dofs());
    postprocessing_data.domain_vm_stress = reduced_smoothing_solution;
  } else
    throw std::runtime_error("Incorrect postprocessing data request.");
}

void ElasticityLin::Run() {
  _SetupSystem();
  _AssembleSystem();
  _AssembleRhs();
  _SolveSystem();
  _Postprocess();
} // end of run

Vector<double> ElasticityLin::GetPostprocessingData(unsigned int &data_flag) {
  Vector<double> postprocessing_vector;
  if (data_flag == e_von_mises_s)
    postprocessing_vector = postprocessing_data.domain_vm_stress;
  else
    throw std::runtime_error("Incorrect postprocessing data request.");
  return postprocessing_vector;
}