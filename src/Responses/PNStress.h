#pragma once

#include <ResponseBase.h>

class PNStress : public ResponseBase {
public:
  PNStress(BVP &, const unsigned int &);
  virtual ~PNStress() = default;

  double GetFunction() override;
  double GetFunctionRef() override;

  // Pure virtual: Must be implemented by Linear/NonLinear classes
  virtual void ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                          double &cell_gradient_density) override = 0;

  virtual void GetOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) override;

protected:
  // Core logic that stays in Base
  void _RunAdjointBVP();
  void _UpdateCorrectionParameters();

  virtual void _ComputeHelpVariables() = 0;
  virtual void _AssembleAdjointRhs() = 0;

  std::map<CellId, double> _el_stress_map;

  double _p_value = 0;
  double _integral_p_norm_vm = 0.0;
  double _stress_limit = 0.0;
  double _corrector = 1.0;
  double _current_max_stress = 0.0;
  double _volume = 0.0;
  double _vol_mult = 0.0;

  LA::MPI::Vector _adjoint_rhs;
  Vector<double> _adjoint_solution;

  // For the output
  Vector<double> _element_stresses;

  unsigned int _opt_iter = 0;
};