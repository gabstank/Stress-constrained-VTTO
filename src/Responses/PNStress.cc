// deal.II files
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

// Project files
#include <DesignField.h>
#include <Logger.h>
#include <PNStress.h>

PNStress::PNStress(BVP &bvp_, const unsigned int &id_) : ResponseBase(bvp_, id_) {
  auto params = parameters();

  _stress_limit = 1.;
  for (const auto &constr : params.optimization.constraints) {
    if (constr.id == this->_id) {
      _stress_limit = constr.constr_value;
      break;
    }
  }
  _p_value = 10;
}

double PNStress::GetFunction() {
  _RunAdjointBVP();
  this->_function_value = _corrector * std::pow(_integral_p_norm_vm / _volume, 1.0 / _p_value);
  auto &logger = Logger::getInstance();
  logger.addData("max_stress", _current_max_stress);
  logger.addData("corrector", _corrector);
  logger.addData("PNfunction", this->_function_value / _corrector);
  ++_opt_iter;
  return this->_function_value;
}

double PNStress::GetFunctionRef() {
  this->_ref_value =
      1.0; // The reference value for the stress constraint can be set to 1.0, as the constraint is
           // defined as the ratio of the current stress to the stress limit. This means that when
           // the function value is equal to 1.0, the current stress is equal to the stress limit,
           // which is the desired outcome for the constraint.
  return this->_ref_value;
}

void PNStress::GetOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {
  std::vector<std::string> stress_name(1, "element_stresses");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      scalar_data_component_interpretation(DataComponentInterpretation::component_is_scalar);
  output.template PushDataName<Vector<double>>(_element_stresses, stress_name,
                                               scalar_data_component_interpretation,
                                               &this->_bvp._dof_handler);
}

void PNStress::_RunAdjointBVP() {
  _UpdateCorrectionParameters();
  _ComputeHelpVariables();
  _adjoint_rhs.reinit(this->_bvp._locally_owned_dofs, mpi_communicator());
  _adjoint_solution.reinit(this->_bvp._dof_handler.n_dofs());
  _AssembleAdjointRhs();
  this->_bvp._SolveLinearSystem(this->_bvp._system_matrix, _adjoint_solution, _adjoint_rhs);

  PROJ_MPI_BARRIER
}

void PNStress::_UpdateCorrectionParameters() {
  // Adaptive stress constraint
  if (_current_max_stress == 0.) {
    _corrector = 1.;
  } else {
    double uncorrected_function_value = std::pow(_integral_p_norm_vm / _volume, 1.0 / _p_value);
    _corrector =
        1. * ((_current_max_stress / _stress_limit) / uncorrected_function_value) + 0. * _corrector;
  }

  // pcout() << "Max stress: " << _current_max_stress;
  // pcout() << " | Function value/c: " << this->_function_value / _corrector;
  // pcout() << " | c: " << _corrector;
  // pcout() << " | stress limit: " << _stress_limit << std::endl;
}
