// Deal.II headers
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

// Project headers
#include <Logger.h>
#include <PDEFilter.h>

PDEFilter::PDEFilter(Mesh &mesh_, std::string type)
    : Filter(mesh_), _type(type), _filter_dof_handler(mesh_.tria), _fe_nothing(FE_Nothing<2>(), 1),
      _fe_filter(FE_Q<2>(1), 1),
      _filter_length_parameter(parameters().optimization.filters.blur_filter_radius /
                               (2.0 * std::sqrt(3.0)))
{
  _fe_collection.push_back(_fe_filter);
  _fe_collection.push_back(_fe_nothing);
}

void PDEFilter::updateContinuationParameters()
{
  // For now, now continuation parameters for PDE filter
}

void PDEFilter::applyFilter(std::map<CellId, double> &density_map)
{
  _in_out_map.clear();
  _in_out_map = density_map;
  auto params = parameters();
  _SetActiveIndex();
  _SetupFilterSystem();
  _AssembleFilterSystem();
  if (params.optimization.filters.pde_boundary_penalization)
    _AssembleFilterSystemBoundary();
  _AssembleFilterRHS();
  _SolveFilterSystem();
  if (_type == "density")
    _InterpolateCellValuesFromNodes(true);
  else
    _InterpolateCellValuesFromNodes(false);
  density_map.clear();
  density_map = _in_out_map;
}

void PDEFilter::applyChainrule(std::map<CellId, double> &sensitivity_map)
{
  auto params = parameters();
  _in_out_map.clear();
  _in_out_map = sensitivity_map;
  _SetActiveIndex();
  _SetupFilterSystem();
  _AssembleFilterSystem();
  if (params.optimization.filters.pde_boundary_penalization)
    _AssembleFilterSystemBoundary();
  _AssembleFilterRHS();
  _SolveFilterSystem();
  _InterpolateCellValuesFromNodes(false);
  sensitivity_map.clear();
  sensitivity_map = _in_out_map;
}

void PDEFilter::_SetActiveIndex()
{
  for (const auto &cell : _filter_dof_handler.active_cell_iterators())
  {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;
    cell->set_active_fe_index(0);
  }
}

void PDEFilter::_SetupFilterSystem()
{
  _filter_dof_handler.distribute_dofs(_fe_collection);

  // Update constraints
  _filter_constraints.clear();
  DoFTools::make_hanging_node_constraints(_filter_dof_handler, _filter_constraints);
  _filter_constraints.close();

  // Update locally owned dofs
  _locally_owned_dofs = _filter_dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(_filter_dof_handler, _locally_relevant_dofs);

  _filter_rhs.reinit(_locally_owned_dofs, mpi_communicator());
  _filtered_field.reinit(_locally_owned_dofs, _locally_relevant_dofs, mpi_communicator());

  // Setup matrix
  DynamicSparsityPattern dsp(_locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(_filter_dof_handler, dsp, _filter_constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp, _locally_owned_dofs, mpi_communicator(),
                                             _locally_relevant_dofs);
  _filter_matrix.reinit(_locally_owned_dofs, _locally_owned_dofs, dsp, mpi_communicator());
}

void PDEFilter::_AssembleFilterSystem()
{
  auto params = parameters();
  _filter_matrix = 0.;

  // Old
  /*QGauss<2> quadrature_formula(_fe_filter.degree + 1);
  FEValues<2> fe_values(_fe_filter, quadrature_formula,
                          update_values | update_gradients | update_quadrature_points |
                              update_JxW_values);

  const unsigned int dofs_per_cell = _fe_filter.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);*/

  // NEW
  const QGauss<2> void_quadrature(_fe_filter.degree + 1);
  const QGauss<2> solid_quadrature(_fe_filter.degree + 1);
  hp::QCollection<2> q_collection;
  q_collection.push_back(solid_quadrature);
  q_collection.push_back(void_quadrature);
  hp::FEValues<2> fe_values_hp(_fe_collection, q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

  unsigned int dofs_per_cell, n_q_points;
  FullMatrix<double> cell_matrix;
  std::vector<types::global_dof_index> local_dof_indices;

  const FEValuesExtractors::Scalar density_fe(0);

  for (const auto &cell : _filter_dof_handler.active_cell_iterators())
  {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;

    // Old
    // fe_values.reinit(cell);

    // New
    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    n_q_points = fe_values.n_quadrature_points;
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);

    cell_matrix = 0.0;

    double r_mod = _filter_length_parameter;

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      double JxW = fe_values.JxW(q);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        Tensor<1, 2> grad_N_i = fe_values[density_fe].gradient(i, q);
        double N_i = fe_values[density_fe].value(i, q);
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          Tensor<1, 2> grad_N_j = fe_values[density_fe].gradient(j, q);
          double N_j = fe_values[density_fe].value(j, q);
          cell_matrix(i, j) += (std::pow(r_mod, 2.0) * grad_N_i * grad_N_j + N_i * N_j) * JxW;

        } // end of j loop
      } // end of i loop
    } // end of quadrature loop

    cell->get_dof_indices(local_dof_indices);

    _filter_constraints.distribute_local_to_global(cell_matrix, local_dof_indices, _filter_matrix);
  } // end of cell loop

  _filter_matrix.compress(VectorOperation::add);
}

void PDEFilter::_AssembleFilterRHS()
{
  _filter_rhs = 0.;

  // Old
  /*QGauss<2> quadrature_formula(_fe_filter.degree + 1);
  FEValues<2> fe_values(_fe_filter, quadrature_formula,
                          update_values | update_gradients | update_quadrature_points |
                              update_JxW_values);

  const unsigned int dofs_per_cell = _fe_filter.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();
  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);*/

  // NEW
  const QGauss<2> void_quadrature(_fe_filter.degree + 1);
  const QGauss<2> solid_quadrature(_fe_filter.degree + 1);
  hp::QCollection<2> q_collection;
  q_collection.push_back(solid_quadrature);
  q_collection.push_back(void_quadrature);
  hp::FEValues<2> fe_values_hp(_fe_collection, q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

  unsigned int dofs_per_cell, n_q_points;
  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  const FEValuesExtractors::Scalar density_fe(0);

  for (const auto &cell : _filter_dof_handler.active_cell_iterators())
  {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;

    // Old
    // fe_values.reinit(cell);
    // cell_rhs = 0.0;

    // New
    fe_values_hp.reinit(cell);
    dofs_per_cell = cell->get_fe().dofs_per_cell;
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    n_q_points = fe_values.n_quadrature_points;

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      double JxW = fe_values.JxW(q);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        double N_i = fe_values[density_fe].value(i, q);
        // Density field constant within a cell, so we can use cell for
        // interpolation
        cell_rhs(i) += _in_out_map.at(cell->id()) * N_i * JxW;
      } // end of i loop
    } // end of quadrature loop

    cell->get_dof_indices(local_dof_indices);

    _filter_constraints.distribute_local_to_global(cell_rhs, local_dof_indices, _filter_rhs);
  } // end of cell loop;
  _filter_rhs.compress(VectorOperation::add);
}

void PDEFilter::_AssembleFilterSystemBoundary()
{

  // Old
  /*QGauss<1> face_quadrature_formula(_fe_filter.degree + 1);
  FEFaceValues<2> fe_face_values(_fe_filter, face_quadrature_formula,
                                   update_values | update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = _fe_filter.dofs_per_cell;
  const unsigned int n_face_q_points = face_quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);*/

  // NEW
  const QGauss<1> face_void_quadrature(_fe_filter.degree + 1);
  const QGauss<1> face_solid_quadrature(_fe_filter.degree + 1);
  hp::QCollection<1> face_q_collection;
  face_q_collection.push_back(face_solid_quadrature);
  face_q_collection.push_back(face_void_quadrature);
  hp::FEFaceValues<2> fe_face_values_hp(_fe_collection, face_q_collection,
                                        update_values | update_quadrature_points |
                                            update_JxW_values);

  unsigned int dofs_per_cell, n_face_q_points;
  FullMatrix<double> cell_matrix;
  std::vector<types::global_dof_index> local_dof_indices;

  const FEValuesExtractors::Scalar density_fe(0);

  for (const auto &cell : _filter_dof_handler.active_cell_iterators())
  {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;

    for (const auto &face : cell->face_iterators())
    {
      if (face->at_boundary() && face->boundary_id() == 0)
      {

        // New
        fe_face_values_hp.reinit(cell, face);
        const FEFaceValues<2> &fe_face_values = fe_face_values_hp.get_present_fe_values();
        n_face_q_points = fe_face_values.n_quadrature_points;
        dofs_per_cell = cell->get_fe().dofs_per_cell;
        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        local_dof_indices.resize(dofs_per_cell);

        // Old
        // fe_face_values.reinit(cell, face);
        // cell_matrix = 0.0;

        for (unsigned int q = 0; q < n_face_q_points; ++q)
        {
          double JxW = fe_face_values.JxW(q);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            double N_i = fe_face_values[density_fe].value(i, q);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              double N_j = fe_face_values[density_fe].value(j, q);

              cell_matrix(i, j) += _filter_length_parameter * N_i * N_j * JxW;
            } // end of j loop
          } // end of i loop
        } // end of quadrature loop

        cell->get_dof_indices(local_dof_indices);

        _filter_constraints.distribute_local_to_global(cell_matrix, local_dof_indices,
                                                       _filter_matrix);
      } // end of if face is at boundary
    } // end of face loop
  } // end of cell loop

  _filter_matrix.compress(VectorOperation::add);

} // end of _AssembleFilterSystemBoundary

void PDEFilter::_SolveFilterSystem()
{
  auto params = parameters();
  SolverControl solver_control(_filter_rhs.size(), params.bvp.solver.tol * _filter_rhs.l2_norm());

  LA::MPI::PreconditionJacobi preconditioner;
  LA::MPI::PreconditionJacobi::AdditionalData additional_data;

  preconditioner.initialize(_filter_matrix, additional_data);

#ifdef USE_PETSC
  LA::SolverCG solver(solver_control, mpi_communicator());
#else
  LA::SolverCG solver(solver_control);
#endif

  LA::MPI::Vector completely_distributed_solution(_filter_rhs);
  completely_distributed_solution = 0.0;
  _filter_constraints.set_zero(completely_distributed_solution);
  solver.solve(_filter_matrix, completely_distributed_solution, _filter_rhs, preconditioner);
  _filter_constraints.distribute(completely_distributed_solution);
  _filtered_field = completely_distributed_solution;
}

void PDEFilter::_InterpolateCellValuesFromNodes(bool set_negative_values_to_zero)
{
  const unsigned int dofs_per_cell = _fe_filter.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : _filter_dof_handler.active_cell_iterators())
  {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;
    // New
    if (cell->active_fe_index() == 1)
      continue;
    _in_out_map[cell->id()] = 0.;

    cell->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      _in_out_map[cell->id()] += 0.25 * _filtered_field[local_dof_indices[i]];
    // Set all negative values to 0
    if (set_negative_values_to_zero && _in_out_map[cell->id()] < 0.)
      _in_out_map[cell->id()] = 0.;
  }
}