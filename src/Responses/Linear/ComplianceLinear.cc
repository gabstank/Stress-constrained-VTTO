// Project files
#include <ComplianceLinear.h>
#include <DesignField.h>

// deal.II Files
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

ComplianceLinear::ComplianceLinear(BVP &bvp_, const unsigned int &id_)
    : ResponseBase(bvp_, id_), _l_mat(std::make_unique<MaterialLinElastic>()) {}

ComplianceLinear::~ComplianceLinear() = default;

double ComplianceLinear::GetFunction() {
  Vector<double> global_external_force(this->_bvp._external_force);
  this->_function_value = this->_bvp._solution * global_external_force;
  return this->_function_value;
}

void ComplianceLinear::ComputeCellDensityGradient(
    const typename DoFHandler<2>::active_cell_iterator &cell, double &cell_gradient_density) {

  auto params = parameters();

  const FEValuesExtractors::Vector u_fe(0);
  hp::FEValues<2> fe_values_hp(this->_bvp._fe_collection, this->_bvp._q_collection,
                               update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
  fe_values_hp.reinit(cell);
  const FEValues<2> &fe_values = fe_values_hp.get_present_fe_values();
  unsigned int n_q_points = fe_values.n_quadrature_points;

  std::vector<Tensor<2, 2>> sol_grads_u(n_q_points);

  fe_values[u_fe].get_function_gradients(this->_bvp._solution, sol_grads_u);

  double cell_density = DesignField::getCurrent().at(cell->id());
  double physical_density_der = DesignField::GetPhysicalDensityDer(cell_density);

  cell_gradient_density = 0.;

  for (unsigned int q = 0; q < n_q_points; ++q) {

    SymmetricTensor<4, 2> elasticity_tensor;
    elasticity_tensor = _l_mat->GetElasticityTensor();
    cell_gradient_density -= symmetrize(sol_grads_u[q]) * physical_density_der * elasticity_tensor *
                             symmetrize(sol_grads_u[q]) * fe_values.JxW(q);

  } // end of loop over q points
}