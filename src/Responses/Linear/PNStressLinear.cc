#include <DesignField.h>
#include <PNStressLinear.h>

PNStressLinear::PNStressLinear(BVP &bvp_, const unsigned int &id_)
    : PNStress(bvp_, id_), _l_mat(std::make_unique<MaterialLinElastic>()) {}

void PNStressLinear::_ComputeHelpVariables() {
  double local_sum_pvm = 0.0;
  double local_vol = 0.0;
  _integral_p_norm_vm = 0.;
  _volume = 0.;
  _current_max_stress = 0.;

  // We are interested in the stress at the centroid of the cell, so we can use a single quadrature
  // point at the centroid for the integration. This is a common approach in topology optimization
  // when evaluating stress-based objective functions, as it provides a good balance between
  // accuracy and computational efficiency. Using more quadrature points can increase the accuracy
  // of the stress evaluation, but it also increases the computational cost. By using a single
  // quadrature point at the centroid, we can capture the average stress state of the cell while
  // keeping the computational cost manageable.
  hp::QCollection<2> q_collection_midpoint;
  q_collection_midpoint.push_back(QMidpoint<2>());

  hp::FEValues<2> fe_values_hp(this->_bvp._fe_collection, q_collection_midpoint,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
  const FEValuesExtractors::Vector u_fe(0);
  auto densities = DesignField::getCurrent();

  // Cell counter for stress output
  int i = 0;
  const unsigned int n_cells = this->_bvp._dof_handler.get_triangulation().n_active_cells();
  _element_stresses.reinit(n_cells);
  _el_stress_map.clear();
  _el_stress_map =
      this->_bvp
          .GetElementStresses(); // Get the (possibly filtered) elemental stresses from the BVP

  for (auto &cell : this->_bvp._dof_handler.active_cell_iterators()) {

    if (!cell->is_locally_owned()) {
      ++i;
      continue;
    }

    fe_values_hp.reinit(cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    unsigned int n_q_points = fe_values.n_quadrature_points;

    for (unsigned int q = 0; q < n_q_points; ++q) {
      double relaxed_stress = _el_stress_map.at(cell->id());
      local_sum_pvm += std::pow(relaxed_stress / _stress_limit, _p_value) * fe_values.JxW(q);
      local_vol += fe_values.JxW(q);

      // We have only one q point so we can increment the cell counter here for the output
      _element_stresses[i] = relaxed_stress;
      ++i;

      if (_el_stress_map.at(cell->id()) > _current_max_stress)
        _current_max_stress = _el_stress_map.at(cell->id());
    } // end of loop over quadrature points
  } // loop over cells

  _integral_p_norm_vm = Utilities::MPI::sum(local_sum_pvm, mpi_communicator());
  _volume = Utilities::MPI::sum(local_vol, mpi_communicator());
  _current_max_stress = Utilities::MPI::max(_current_max_stress, mpi_communicator());
}

void PNStressLinear::_AssembleAdjointRhs() {
  this->_adjoint_rhs = 0.0;

  hp::QCollection<2> q_collection_midpoint;
  q_collection_midpoint.push_back(QMidpoint<2>());

  hp::FEValues<2> fe_values_hp(this->_bvp._fe_collection, q_collection_midpoint,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
  const FEValuesExtractors::Vector u_fe(0);

  SymmetricTensor<4, 2> elasticity_tensor = _l_mat->GetElasticityTensor();
  double outer_part = std::pow(_integral_p_norm_vm / _volume, (1. / _p_value) - 1.) / _volume;
  auto densities = DesignField::getCurrent();

  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  // actual run to assemble rhs vector
  for (const auto &cell : this->_bvp._dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    double cell_density = densities.at(cell->id());
    double cell_stress_relaxation = DesignField::StressRelaxation(cell_density);

    fe_values_hp.reinit(cell);
    const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
    unsigned int n_q_points = fe_values.n_quadrature_points;
    unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);

    std::vector<SymmetricTensor<2, 2>> solution_sym_grads_u(fe_values.n_quadrature_points);
    fe_values[u_fe].get_function_symmetric_gradients(this->_bvp._solution, solution_sym_grads_u);

    for (unsigned int q = 0; q < n_q_points; ++q) {
      double vm_stress = _l_mat->GetVonMisesStress(solution_sym_grads_u[q]);
      SymmetricTensor<2, 2> cauchy_stress = _l_mat->GetCauchyStress(solution_sym_grads_u[q]);
      double relaxed_stress = cell_stress_relaxation * vm_stress;
      double stress_measure_q = relaxed_stress / _stress_limit;

      double stress_common_der_part = std::pow(stress_measure_q, _p_value - 1.) / _stress_limit *
                                      0.5 / vm_stress * cell_stress_relaxation;

      if (std::isnan(stress_common_der_part) || std::isinf(stress_common_der_part) ||
          stress_common_der_part == 0.)
        continue;

      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        SymmetricTensor<2, 2> Grad_N_i = fe_values[u_fe].symmetric_gradient(i, q);
        SymmetricTensor<2, 2> der_cauchy_stress_prt1 = Grad_N_i * elasticity_tensor;

        cell_rhs(i) -= -1.0 * (cauchy_stress * StandardTensors::I()) *
                       (der_cauchy_stress_prt1 * StandardTensors::I()) * fe_values.JxW(q) *
                       stress_common_der_part;

        cell_rhs(i) -= 3.0 * der_cauchy_stress_prt1 * cauchy_stress * fe_values.JxW(q) *
                       stress_common_der_part;

      } // end of loop over i

    } // end of loop over q points

    cell_rhs *= outer_part;

    cell->get_dof_indices(local_dof_indices);
    this->_bvp._constraints.distribute_local_to_global(cell_rhs, local_dof_indices, _adjoint_rhs);
  }
  _adjoint_rhs.compress(VectorOperation::add);
}

void PNStressLinear::ComputeCellDensityGradient(
    const typename DoFHandler<2>::active_cell_iterator &cell, double &cell_gradient_density) {
  const FEValuesExtractors::Vector u_fe(0);

  hp::QCollection<2> q_collection_midpoint;
  q_collection_midpoint.push_back(QMidpoint<2>());

  hp::FEValues<2> fe_values_hp(this->_bvp._fe_collection, q_collection_midpoint,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
  fe_values_hp.reinit(cell);
  const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
  unsigned int n_q_points = fe_values.n_quadrature_points;

  std::vector<SymmetricTensor<2, 2>> sol_sym_grads_u(n_q_points);
  std::vector<SymmetricTensor<2, 2>> adj_sol_sym_grads_u(n_q_points);

  fe_values[u_fe].get_function_symmetric_gradients(this->_bvp._solution, sol_sym_grads_u);
  fe_values[u_fe].get_function_symmetric_gradients(_adjoint_solution, adj_sol_sym_grads_u);

  double cell_density = DesignField::getCurrent().at(cell->id());
  double physical_density_der = DesignField::GetPhysicalDensityDer(cell_density);
  double cell_stress_relaxation = DesignField::StressRelaxation(cell_density);
  double cell_stress_relaxation_der = DesignField::StressRelaxationDer(cell_density);

  double outer_part = std::pow(_integral_p_norm_vm / _volume, (1. / _p_value) - 1.) / _volume;
  cell_gradient_density = 0.0;
  SymmetricTensor<4, 2> elasticity_tensor = _l_mat->GetElasticityTensor();

  for (unsigned int q = 0; q < n_q_points; ++q) {
    double vm_stress = _l_mat->GetVonMisesStress(sol_sym_grads_u[q]);
    double relaxed_stress = cell_stress_relaxation * vm_stress;
    double stress_measure_q = relaxed_stress / _stress_limit;
    double stress_relaxation_der_part = std::pow(stress_measure_q, _p_value - 1.) / _stress_limit *
                                        cell_stress_relaxation_der * vm_stress * fe_values.JxW(q);

    double adjoint_part = adj_sol_sym_grads_u[q] * physical_density_der * elasticity_tensor *
                          sol_sym_grads_u[q] * fe_values.JxW(q);

    cell_gradient_density += stress_relaxation_der_part * outer_part;
    cell_gradient_density += adjoint_part;
  }
  cell_gradient_density *= _corrector;
}