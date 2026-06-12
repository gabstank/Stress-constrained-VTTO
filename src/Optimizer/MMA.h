#pragma once

// Project headers
#include <Optimizer.h>

// deal.II headers

using namespace dealii;

class MMA : public Optimizer {

  double _fscale = 0.0;
  const double _movlim;
  // "speed-control" for the asymptotes
  double _asyminit, _asymdec, _asyminc;

  std::map<unsigned int, double> _response_values;
  std::map<unsigned int, std::map<CellId, double>> _gradients;

  std::map<CellId, double> _pseudo_densities;
  std::map<CellId, double> _old_pseudo_densities;
  std::map<CellId, double> _old2_pseudo_densities;
  std::map<CellId, double> _pseudo_densities_min;
  std::map<CellId, double> _pseudo_densities_max;
  std::map<CellId, double> _low_asym, _upp_asym; // Asymptotes
  std::map<CellId, double> _alpha, _beta;        // Asymptotes bounds
  std::map<CellId, double> _p0, _q0;
  std::map<unsigned int, std::map<CellId, double>> _pij, _qij;
  // Local: subproblem constant terms, dual gradient, dual hessian
  std::map<unsigned int, double> _b, _grad;
  std::map<unsigned int, std::map<unsigned int, double>> _hess;
  // Local vectors: penalty numbers for subproblem
  std::map<unsigned int, double> _a, _c, _d, _lam, _mu, _s, _y;
  double _z = 0.0;
  double _epsi = 1.0;

  unsigned int _iteration = 0;

  void _SetOuterMoveLimit();
  void _GenSub();
  void _SolveDIP();
  void _XYZofLAMBDA();
  void _DualGrad();
  void _DualHess();
  void _Solve();
  void _DualLineSearch();
  double _DualResidual();

public:
  MMA();
  ~MMA() override = default;

  void DesignUpdate(std::map<CellId, double> &pseudo_densities,
                    const std::map<unsigned int, std::map<CellId, double>> &gradients,
                    const std::map<unsigned int, double> &response_values) override;

  std::map<std::string, std::map<CellId, double>> GetDataRelevantForTransfer() override;
  void SetDataAfterTransfer(const std::map<CellId, double> &data, const std::string &name) override;
};