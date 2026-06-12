// C++ headers

// Deal.II headers
#include <deal.II/fe/fe_series.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/smoothness_estimator.h>

// Project headers
#include <BVP.h>

BVP::BVP(Mesh &mesh_)
    : _smoothing_fe(FE_Q<2>(parameters().bvp.poly_degree)), _smoothing_dof_handler(mesh_.tria),
      _tria(mesh_.tria), _dof_handler(mesh_.tria), _stress_filter(mesh_.tria) {

  _fe_collection.push_back(FESystem<2>(FE_Q<2>(parameters().bvp.poly_degree), 2));
  _q_collection.push_back(QGauss<2>(parameters().bvp.poly_degree + 1));
  _face_q_collection.push_back(QGauss<1>(parameters().bvp.poly_degree + 1));

  _smoothing_fe_collection.push_back(FE_Q<2>(parameters().bvp.poly_degree));
  _smoothing_q_collection.push_back(QGauss<2>(parameters().bvp.poly_degree + 1));
}

void BVP::_SolveLinearSystem(LA::MPI::SparseMatrix &system_matrix, Vector<double> &solution,
                             LA::MPI::Vector &system_rhs) {

  TimerOutput::Scope t(compute_timer(), "BVP::_SolveLinearSystem");

  auto params = parameters();

  LA::MPI::Vector completely_distributed_solution(system_rhs);

  completely_distributed_solution = 0.0;
  _constraints.set_zero(completely_distributed_solution);

  SolverControl solver_control(system_rhs.size(), params.bvp.solver.tol * system_rhs.l2_norm());

  LA::MPI::PreconditionJacobi preconditioner;
  LA::MPI::PreconditionJacobi::AdditionalData additional_data;

  preconditioner.initialize(system_matrix, additional_data);

  if (params.bvp.solver.name == "CG") {
#ifdef USE_PETSC
    LA::SolverCG solver(solver_control, mpi_communicator());
#else
    LA::SolverCG solver(solver_control);
#endif
    solver.solve(system_matrix, completely_distributed_solution, system_rhs, preconditioner);

  } else if (params.bvp.solver.name == "Direct") {

#ifdef USE_PETSC
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator());
    solver.initialize(preconditioner);
    solver.set_symmetric_mode(true);
#else
    TrilinosWrappers::SolverDirect::AdditionalData additional_data;
    additional_data.solver_type = "Amesos_Superludist";
    TrilinosWrappers::SolverDirect solver(solver_control, additional_data);
#endif
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
  } else {
    throw std::runtime_error("Wrong linear solver type.");
  }

  _constraints.distribute(completely_distributed_solution);
  solution = completely_distributed_solution;
}

void BVP::_SetupSmoothingSystem() {

  auto params = parameters();

  _smoothing_dof_handler.distribute_dofs(_smoothing_fe_collection);

  _locally_owned_smooth_dofs = _smoothing_dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(_smoothing_dof_handler, _locally_relevant_smooth_dofs);

  _smoothing_solution.reinit(_locally_owned_smooth_dofs, _locally_relevant_smooth_dofs,
                             mpi_communicator());
  _smoothing_rhs.reinit(_locally_owned_smooth_dofs, mpi_communicator());

  _smoothing_constraints.clear();
  _smoothing_constraints.reinit(_locally_relevant_smooth_dofs);
  DoFTools::make_hanging_node_constraints(_smoothing_dof_handler, _smoothing_constraints);
  _smoothing_constraints.close();

  DynamicSparsityPattern dsp(_locally_relevant_smooth_dofs);
  DoFTools::make_sparsity_pattern(_smoothing_dof_handler, dsp, _smoothing_constraints, true);
  SparsityTools::distribute_sparsity_pattern(dsp, _locally_owned_smooth_dofs, mpi_communicator(),
                                             _locally_relevant_smooth_dofs);

  _smoothing_matrix.reinit(_locally_owned_smooth_dofs, _locally_owned_smooth_dofs, dsp,
                           mpi_communicator());
}

void BVP::AssembleSmoothingMatrix() {

  _smoothing_matrix = 0;

  // QGauss<2> quadrature_formula(_smoothing_fe.degree + 1);
  hp::FEValues<2> fe_values_hp(_smoothing_fe_collection, _smoothing_q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

  double n_i, n_j;

  FullMatrix<double> cell_matrix;
  std::vector<unsigned int> local_dof_indices;
  unsigned int dofs_per_cell;

  for (auto &cell : _smoothing_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);

    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
      // phi_phi assembly
      for (unsigned int i = 0; i < dofs_per_cell; i++) {
        n_i = fe_values.shape_value(i, q);
        for (unsigned int j = 0; j < dofs_per_cell; j++) {
          n_j = fe_values.shape_value(j, q);
          cell_matrix(i, j) += n_i * n_j * fe_values.JxW(q);
        } // end of loop over j
      } // end of loop over i

    } // end of loop over q points

    // add to system smoothing matrix
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
    _smoothing_constraints.distribute_local_to_global(cell_matrix, local_dof_indices,
                                                      _smoothing_matrix);
  }
  _smoothing_matrix.compress(VectorOperation::add);
}

void BVP::SolveSmoothingSystem() {

  auto params = parameters();

  SolverControl solver_control(_smoothing_dof_handler.n_dofs(),
                               params.bvp.solver.tol * _smoothing_rhs.l2_norm());

  LA::MPI::PreconditionJacobi preconditioner;
  LA::MPI::PreconditionJacobi::AdditionalData additional_data;

  preconditioner.initialize(_smoothing_matrix, additional_data);

  LA::MPI::Vector completely_distributed_solution(_smoothing_rhs);
  completely_distributed_solution = 0.0;

  _smoothing_constraints.set_zero(completely_distributed_solution);

  if (params.bvp.solver.name == "CG") {
#ifdef USE_PETSC
    LA::SolverCG solver(solver_control, mpi_communicator());
#else
    LA::SolverCG solver(solver_control);
#endif
    solver.solve(_smoothing_matrix, completely_distributed_solution, _smoothing_rhs,
                 preconditioner);
  } else {
#ifdef USE_PETSC
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator());
    solver.initialize(preconditioner);
    solver.set_symmetric_mode(true);
#else
    TrilinosWrappers::SolverDirect::AdditionalData additional_data;
    additional_data.solver_type = "Amesos_Superludist";
    TrilinosWrappers::SolverDirect solver(solver_control, additional_data);
#endif
    solver.solve(_smoothing_matrix, completely_distributed_solution, _smoothing_rhs);
  }

  _smoothing_constraints.distribute(completely_distributed_solution);

  _smoothing_solution = completely_distributed_solution;
}

void BVP::GetAllOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {

  auto di = DataComponentInterpretationTypes::vec_interpretation;
  auto di_scalar = DataComponentInterpretationTypes::scalar_interpretation;

  std::vector<std::string> sol_names(2, "sol_u"), rhs_names(2, "rhs_u"), mat_names(1, "material"),
      el_stress_raw_names(1, "el_stress_raw"), el_stress_smoothed_names(1, "el_stress_smoothed");

  Vector<double> materialid(this->_tria.n_active_cells());
  Vector<double> el_stress_raw(this->_tria.n_active_cells());
  Vector<double> el_stress_smoothed(this->_tria.n_active_cells());
  int i = 0;
  for (auto cell : this->_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) {
      ++i;
      continue;
    }
    materialid[i] = cell->material_id();
    el_stress_raw[i] = _el_stress_map_raw.at(cell->id());
    el_stress_smoothed[i] = _el_stress_map_smoothed.at(cell->id());
    ++i;
  }

  output.template PushDataName<LA::MPI::Vector>(_system_rhs, rhs_names, di, &_dof_handler);
  output.template PushDataName<Vector<double>>(_solution, sol_names, di, &_dof_handler);
  output.template PushDataName<Vector<double>>(materialid, mat_names, di_scalar,
                                               &this->_smoothing_dof_handler);
  output.template PushDataName<Vector<double>>(el_stress_raw, el_stress_raw_names, di_scalar,
                                               &this->_smoothing_dof_handler);
  output.template PushDataName<Vector<double>>(el_stress_smoothed, el_stress_smoothed_names,
                                               di_scalar, &this->_smoothing_dof_handler);
} // end of function

Vector<double> BVP::GetRefinementData(const std::string &ref_type) {
  (void)ref_type;
  Vector<double> vec;
  throw std::runtime_error("Wrong refinement type.");
  return vec;
} // end of function