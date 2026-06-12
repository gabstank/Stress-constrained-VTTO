// Project files
#include <MMA.h>
#include <Parameter.h>

MMA::MMA()
    : _movlim(parameters().optimization.mma.move_limit),
      _asyminit(parameters().optimization.mma.asym_move_init),
      _asymdec(parameters().optimization.mma.asym_move_dec),
      _asyminc(parameters().optimization.mma.asym_move_incr) {
  for (unsigned int i = 0; i < parameters().optimization.constraints.size(); ++i) {
    unsigned int cons_id = parameters().optimization.constraints[i].id;
    _a[cons_id] = 0.0;
    _c[cons_id] = 1000.0;
    _d[cons_id] = 0.0;
  }
}

void MMA::DesignUpdate(std::map<CellId, double> &pseudo_densities,
                       const std::map<unsigned int, std::map<CellId, double>> &gradients,
                       const std::map<unsigned int, double> &response_values) {
  TimerOutput::Scope t(compute_timer(), "MMA::DesignUpdate");

  _pseudo_densities.clear();
  _pseudo_densities = pseudo_densities;
  _gradients.clear();
  _gradients = gradients;
  _response_values.clear();
  _response_values = response_values;

  _SetOuterMoveLimit();

  // Generate the subproblem
  _GenSub();

  // Update old pseudo density
  _old2_pseudo_densities.clear();
  _old2_pseudo_densities = _old_pseudo_densities;
  _old_pseudo_densities.clear();
  _old_pseudo_densities = _pseudo_densities;

  // Solve the dual with an interior point method
  _SolveDIP();

  pseudo_densities.clear();
  pseudo_densities = _pseudo_densities;

  ++_iteration;

  // Setting pseudo density of non design cells to one.
  /*for(auto &non_design_cell : _non_design_cells){
    if(_pseudo_densities.find(non_design_cell) != _pseudo_densities.end() )
      _pseudo_densities[non_design_cell] = 1.0;
  }*/

} // end of _MMAUpdate()

std::map<std::string, std::map<CellId, double>> MMA::GetDataRelevantForTransfer() {
  std::map<std::string, std::map<CellId, double>> data_for_transfer;
  data_for_transfer[MMATransferDataNames::OldPseudoDensity] = _old_pseudo_densities;
  data_for_transfer[MMATransferDataNames::Old2PseudoDensity] = _old2_pseudo_densities;
  data_for_transfer[MMATransferDataNames::MinDensity] = _pseudo_densities_min;
  data_for_transfer[MMATransferDataNames::MaxDensity] = _pseudo_densities_max;
  data_for_transfer[MMATransferDataNames::LowerAsymptotes] = _low_asym;
  data_for_transfer[MMATransferDataNames::UpperAsymptotes] = _upp_asym;
  return data_for_transfer;
}

void MMA::SetDataAfterTransfer(const std::map<CellId, double> &data, const std::string &name) {
  if (name == MMATransferDataNames::OldPseudoDensity) {
    _old_pseudo_densities.clear();
    _old_pseudo_densities = data;
  } else if (name == MMATransferDataNames::Old2PseudoDensity) {
    _old2_pseudo_densities.clear();
    _old2_pseudo_densities = data;
  } else if (name == MMATransferDataNames::MinDensity) {
    _pseudo_densities_min.clear();
    _pseudo_densities_min = data;
  } else if (name == MMATransferDataNames::MaxDensity) {
    _pseudo_densities_max.clear();
    _pseudo_densities_max = data;
  } else if (name == MMATransferDataNames::LowerAsymptotes) {
    _low_asym.clear();
    _low_asym = data;
  } else if (name == MMATransferDataNames::UpperAsymptotes) {
    _upp_asym.clear();
    _upp_asym = data;
  } else {
    throw std::runtime_error("Invalid data write from MMA class");
  }
}

void MMA::_SetOuterMoveLimit() {
  TimerOutput::Scope t(compute_timer(), "MMA::_SetOuterMoveLimit");

  const double max_density = 1.0, min_density = 0.0;
  for (const auto &[cell_id, value] : _pseudo_densities) {
    _pseudo_densities_max[cell_id] = std::min(max_density, _pseudo_densities[cell_id] + _movlim);
    _pseudo_densities_min[cell_id] = std::max(min_density, _pseudo_densities[cell_id] - _movlim);
  }
} // end of _SetOuterMoveLimit()

void MMA::_GenSub() {

  auto params = parameters();

  if (_iteration < 3) {
    for (const auto &[cell_id, value] : _pseudo_densities) {
      _low_asym[cell_id] =
          _pseudo_densities[cell_id] -
          _asyminit * (_pseudo_densities_max[cell_id] - _pseudo_densities_min[cell_id]);

      _upp_asym[cell_id] =
          _pseudo_densities[cell_id] +
          _asyminit * (_pseudo_densities_max[cell_id] - _pseudo_densities_min[cell_id]);

    } // end of loop over cells
  } // end of if _opt_iteration < 3
  else {

    for (const auto &[cell_id, value] : _pseudo_densities) {
      double gamma, helpvar;

      helpvar = (_pseudo_densities[cell_id] - _old_pseudo_densities[cell_id]) *
                (_old_pseudo_densities[cell_id] - _old2_pseudo_densities[cell_id]);
      if (helpvar < 0.0)
        gamma = _asymdec;
      else if (helpvar > 0.0)
        gamma = _asyminc;
      else
        gamma = 1.0;
      //////////////
      _low_asym[cell_id] = _pseudo_densities[cell_id] -
                           gamma * (_old_pseudo_densities[cell_id] - _low_asym[cell_id]);
      _upp_asym[cell_id] = _pseudo_densities[cell_id] +
                           gamma * (_upp_asym[cell_id] - _old_pseudo_densities[cell_id]);

      double xmi, xma;
      xmi = std::max(1.0e-5, _pseudo_densities_max[cell_id] - _pseudo_densities_min[cell_id]);
      if (params.optimization.mma.robust_asymptotes_type == 0) {
        _low_asym[cell_id] = std::max(_low_asym[cell_id], _pseudo_densities[cell_id] - 10.0 * xmi);
        _low_asym[cell_id] = std::min(_low_asym[cell_id], _pseudo_densities[cell_id] - 0.01 * xmi);
        _upp_asym[cell_id] = std::max(_upp_asym[cell_id], _pseudo_densities[cell_id] + 0.01 * xmi);
        _upp_asym[cell_id] = std::min(_upp_asym[cell_id], _pseudo_densities[cell_id] + 10.0 * xmi);
      } else if (params.optimization.mma.robust_asymptotes_type == 1) {
        _low_asym[cell_id] = std::max(_low_asym[cell_id], _pseudo_densities[cell_id] - 100.0 * xmi);
        _low_asym[cell_id] =
            std::min(_low_asym[cell_id], _pseudo_densities[cell_id] - 1.0e-4 * xmi);
        _upp_asym[cell_id] =
            std::max(_upp_asym[cell_id], _pseudo_densities[cell_id] + 1.0e-4 * xmi);
        _upp_asym[cell_id] = std::min(_upp_asym[cell_id], _pseudo_densities[cell_id] + 100.0 * xmi);
        xmi = _pseudo_densities_min[cell_id] - 1.0e-5;
        xma = _pseudo_densities_max[cell_id] + 1.0e-5;
        if (_pseudo_densities[cell_id] < xmi) {
          _low_asym[cell_id] =
              _pseudo_densities[cell_id] - (xma - _pseudo_densities[cell_id]) / 0.9;
          _upp_asym[cell_id] =
              _pseudo_densities[cell_id] + (xma - _pseudo_densities[cell_id]) / 0.9;
        }
        if (_pseudo_densities[cell_id] > xma) {
          _low_asym[cell_id] =
              _pseudo_densities[cell_id] - (_pseudo_densities[cell_id] - xmi) / 0.9;
          _upp_asym[cell_id] =
              _pseudo_densities[cell_id] + (_pseudo_densities[cell_id] - xmi) / 0.9;
        }
      } else {
        std::string err_msg = std::string(__FILE__) + ":" + std::to_string(__LINE__) +
                              " Incorrect Robust asymptotes type";
        throw std::runtime_error(err_msg);
      }

    } // end of loop over cells
  } // end of if _opt_iteration > 2

  _alpha.clear();
  _beta.clear();
  _p0.clear();
  _q0.clear();
  _pij.clear();
  _qij.clear();

  // double raa0 = 0.5*1e-6;
  const double albefa = 0.1, raa0 = 0.00001, move = 0.5, xmamieps = 1.0e-5;
  unsigned int obj_idx = params.optimization.objective.id;
  // Set bounds and the coefficients for the approximation
  for (const auto &[cell_id, value] : _pseudo_densities) {
    // TODO: PETSC TOP OPT
    /*_alpha[cell_id] = std::max(_pseudo_densities_min[cell_id],
                               0.9 * _low_asym[cell_id] + albefa * _pseudo_densities[cell_id]);
    _beta[cell_id] = std::min(_pseudo_densities_max[cell_id],
                              0.9 * _upp_asym[cell_id] - albefa * _pseudo_densities[cell_id]);*/
    // Compute bounds alpha and beta (Svanberg 2002)
    _alpha[cell_id] =
        std::max(_pseudo_densities_min[cell_id],
                 _low_asym[cell_id] + albefa * (_pseudo_densities[cell_id] - _low_asym[cell_id]));
    _alpha[cell_id] = std::max(_alpha[cell_id], _pseudo_densities[cell_id] -
                                                    move * (_pseudo_densities_max[cell_id] -
                                                            _pseudo_densities_min[cell_id]));
    _alpha[cell_id] = std::min(_alpha[cell_id], _pseudo_densities_max[cell_id]);

    _beta[cell_id] =
        std::min(_pseudo_densities_max[cell_id],
                 _upp_asym[cell_id] - albefa * (_upp_asym[cell_id] - _pseudo_densities[cell_id]));
    _beta[cell_id] = std::min(_beta[cell_id],
                              _pseudo_densities[cell_id] + move * (_pseudo_densities_max[cell_id] -
                                                                   _pseudo_densities_min[cell_id]));
    _beta[cell_id] = std::max(_beta[cell_id], _pseudo_densities_min[cell_id]);

    // Objective function
    double dfdxp = std::max(0.0, _gradients[obj_idx][cell_id]);
    double dfdxm = std::max(0.0, -1.0 * _gradients[obj_idx][cell_id]);

    // (Svanberg 2002)
    double xmamiinv =
        1.0 / std::max(xmamieps, _pseudo_densities_max[cell_id] - _pseudo_densities_min[cell_id]);
    double pq = 0.001 * std::abs(_gradients[obj_idx][cell_id]) + raa0 * xmamiinv;
    _p0[cell_id] = std::pow(_upp_asym[cell_id] - _pseudo_densities[cell_id], 2.0) * (dfdxp + pq);
    _q0[cell_id] = std::pow(_pseudo_densities[cell_id] - _low_asym[cell_id], 2.0) * (dfdxm + pq);

    // Constraints
    for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
      unsigned int cons_id = params.optimization.constraints[i].id;
      double dgdxp = std::max(0.0, _gradients[cons_id][cell_id]);
      double dgdxm = std::max(0.0, -1.0 * _gradients[cons_id][cell_id]);

      double pq = 0.001 * std::abs(_gradients[cons_id][cell_id]) + raa0 * xmamiinv;
      if (params.optimization.mma.constraint_modification) {
        _pij[cons_id][cell_id] =
            std::pow(_upp_asym[cell_id] - _pseudo_densities[cell_id], 2.0) * (dgdxp + pq);
        _qij[cons_id][cell_id] =
            std::pow(_pseudo_densities[cell_id] - _low_asym[cell_id], 2.0) * (dgdxm + pq);
      } else {
        _pij[cons_id][cell_id] =
            std::pow(_upp_asym[cell_id] - _pseudo_densities[cell_id], 2.0) * dgdxp;
        _qij[cons_id][cell_id] =
            std::pow(_pseudo_densities[cell_id] - _low_asym[cell_id], 2.0) * dgdxm;
      }
    }
  }

  // The constant for the constraints
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _b[cons_id] = 0.0;
    for (const auto &[cell_id, value] : _pseudo_densities) {
      _b[cons_id] += _pij[cons_id][cell_id] / (_upp_asym[cell_id] - _pseudo_densities[cell_id]) +
                     _qij[cons_id][cell_id] / (_pseudo_densities[cell_id] - _low_asym[cell_id]);
    }
  }
  // All reduce
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _b[cons_id] = Utilities::MPI::sum(_b[cons_id], mpi_communicator());
  } // end of loop over constraints
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _b[cons_id] += -_response_values[cons_id];
  }

} // end of _GenSub()

void MMA::_SolveDIP() {

  auto params = parameters();

  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _lam[cons_id] = _c[cons_id] / 2.0;
    _mu[cons_id] = 1.0;
  }
  unsigned int n_cells = Utilities::MPI::sum(_pseudo_densities.size(), mpi_communicator());
  double tol = 1.0e-9 * std::sqrt(params.optimization.constraints.size() + n_cells);
  _epsi = 1.0;
  double err = 1.0;
  unsigned int loop = 0;

  while (_epsi > tol) {

    loop = 0;
    while (err > 0.9 * _epsi && loop < 100) {
      loop++;

      _XYZofLAMBDA();
      _DualGrad();

      for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
        unsigned int cons_id = params.optimization.constraints[i].id;
        _grad[cons_id] = -1.0 * _grad[cons_id] - _epsi / _lam[cons_id];
      }
      _DualHess();

      _Solve();

      const unsigned int m = params.optimization.constraints.size() *
                             1000; // this varible has no significance this is just for indexing
      for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
        unsigned int cons_id = params.optimization.constraints[i].id;
        _s[cons_id] = _grad[cons_id];
        _s[m + cons_id] =
            -_mu[cons_id] + _epsi / _lam[cons_id] - _s[cons_id] * _mu[cons_id] / _lam[cons_id];
      }

      _DualLineSearch();
      _XYZofLAMBDA();
      err = _DualResidual();

    } // end of while
    _epsi = _epsi * 0.1;
  } // end of while

} // end of _SolveDIP()

void MMA::_XYZofLAMBDA() {

  auto params = parameters();

  double lamai = 0.0;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    if (_lam[cons_id] < 0.0) {
      _lam[cons_id] = 0;
    }
    _y[cons_id] =
        std::max(0.0,
                 _lam[cons_id] - _c[cons_id]); // Note y=(lam-c)/d - however d is fixed at one !!
    lamai += _lam[cons_id] * _a[cons_id];
  } // end of loop over constraints
  _z = std::max(0.0, 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0
  double pjlam, qjlam;

  for (const auto &[cell_id, value] : _pseudo_densities) {
    pjlam = _p0[cell_id];
    qjlam = _q0[cell_id];
    for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
      unsigned int cons_id = params.optimization.constraints[i].id;
      pjlam += _pij[cons_id][cell_id] * _lam[cons_id];
      qjlam += _qij[cons_id][cell_id] * _lam[cons_id];
    } // end of loop over constraints
    _pseudo_densities[cell_id] =
        (std::sqrt(pjlam) * _low_asym[cell_id] + std::sqrt(qjlam) * _upp_asym[cell_id]) /
        (std::sqrt(pjlam) + std::sqrt(qjlam));

    if (_pseudo_densities[cell_id] < _alpha[cell_id]) {
      _pseudo_densities[cell_id] = _alpha[cell_id];
    }
    if (_pseudo_densities[cell_id] > _beta[cell_id]) {
      _pseudo_densities[cell_id] = _beta[cell_id];
    }

  } // end of loop over cells

} // end of _XYZofLAMBDA()

void MMA::_DualGrad() {
  auto params = parameters();
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _grad[cons_id] = 0.0;

    for (const auto &[cell_id, value] : _pseudo_densities) {
      _grad[cons_id] += _pij[cons_id][cell_id] / (_upp_asym[cell_id] - _pseudo_densities[cell_id]) +
                        _qij[cons_id][cell_id] / (_pseudo_densities[cell_id] - _low_asym[cell_id]);
    } // end of loop over cells
  } // end of loop over constraints

  // all reduce
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _grad[cons_id] = Utilities::MPI::sum(_grad[cons_id], mpi_communicator());
    // pcout() << __LINE__ << " | _grad[cons_id] = " << _grad[cons_id] << std::endl;
  } // end of loop over constraints

  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _grad[cons_id] += -_b[cons_id] - _a[cons_id] * _z - _y[cons_id];
  } // end of loop over constraints

} // end of _DualGrad()

void MMA::_DualHess() {

  auto params = parameters();

  std::map<CellId, double> df2;
  std::map<unsigned int, std::map<CellId, double>> PQ;

  for (const auto &[cell_id, value] : _pseudo_densities) {
    double pjlam = _p0[cell_id];
    double qjlam = _q0[cell_id];
    for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
      unsigned int cons_id = params.optimization.constraints[i].id;
      pjlam += _pij[cons_id][cell_id] * _lam[cons_id];
      qjlam += _qij[cons_id][cell_id] * _lam[cons_id];
      PQ[cons_id][cell_id] =
          _pij[cons_id][cell_id] / std::pow(_upp_asym[cell_id] - _pseudo_densities[cell_id], 2.0) -
          _qij[cons_id][cell_id] / std::pow(_pseudo_densities[cell_id] - _low_asym[cell_id], 2.0);
    }
    df2[cell_id] =
        -1.0 / (2.0 * pjlam / std::pow(_upp_asym[cell_id] - _pseudo_densities[cell_id], 3.0) +
                2.0 * qjlam / std::pow(_pseudo_densities[cell_id] - _low_asym[cell_id], 3.0));
    double xp = (std::sqrt(pjlam) * _low_asym[cell_id] + std::sqrt(qjlam) * _upp_asym[cell_id]) /
                (std::sqrt(pjlam) + std::sqrt(qjlam));
    if (xp < _alpha[cell_id]) {
      df2[cell_id] = 0.0;
    }
    if (xp > _beta[cell_id]) {
      df2[cell_id] = 0.0;
    }
  } // end of loop over cells

  std::map<unsigned int, std::map<CellId, double>> tmp;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    for (const auto &[cell_id, value] : _pseudo_densities) {
      tmp[cons_id][cell_id] = 0.0;
      tmp[cons_id][cell_id] += PQ[cons_id][cell_id] * df2[cell_id];
    } // end of loop over cells
  } // end of loop over constraints

  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id_i = params.optimization.constraints[i].id;
    for (unsigned int j = 0; j < params.optimization.constraints.size(); ++j) {
      unsigned int cons_id_j = params.optimization.constraints[j].id;
      _hess[cons_id_i][cons_id_j] = 0.0;
      for (const auto &[cell_id, value] : _pseudo_densities) {
        _hess[cons_id_i][cons_id_j] += tmp[cons_id_i][cell_id] * PQ[cons_id_j][cell_id];
      } // end of loop over cells
    } // i loop over constraints
  } // j loop over constraints

  // all reduce
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id_i = params.optimization.constraints[i].id;
    for (unsigned int j = 0; j < params.optimization.constraints.size(); ++j) {
      unsigned int cons_id_j = params.optimization.constraints[j].id;
      _hess[cons_id_i][cons_id_j] =
          Utilities::MPI::sum(_hess[cons_id_i][cons_id_j], mpi_communicator());
    }
  }

  double lamai = 0.0;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    if (_lam[cons_id] < 0.0) {
      _lam[cons_id] = 0.0;
    }
    lamai += _lam[cons_id] * _a[cons_id];
    if (_lam[cons_id] > _c[cons_id]) {
      _hess[cons_id][cons_id] += -1.0;
    }
    _hess[cons_id][cons_id] += -_mu[cons_id] / _lam[cons_id];
  } // end of loop over constraints

  if (lamai > 0.0) {
    for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
      unsigned int cons_id_i = params.optimization.constraints[i].id;
      for (unsigned int j = 0; j < params.optimization.constraints.size(); ++j) {
        unsigned int cons_id_j = params.optimization.constraints[j].id;
        _hess[cons_id_i][cons_id_j] += -10.0 * _a[cons_id_i] * _a[cons_id_j];
      } // i loop over constraints
    } // j loop over constraints
  } // lamai > 0.0

  double hess_trace = 0.0;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    hess_trace += _hess[cons_id][cons_id];
  }
  double hess_corr = 1e-4 * hess_trace / params.optimization.constraints.size();
  if (-1.0 * hess_corr < 1.0e-7) {
    hess_corr = -1.0e-7;
  }
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _hess[cons_id][cons_id] += hess_corr;
  }
} // end of _DualHess()

void MMA::_Solve() {

  auto params = parameters();

  FullMatrix<double> hessian(params.optimization.constraints.size());
  Vector<double> grad(params.optimization.constraints.size());

  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id_i = params.optimization.constraints[i].id;
    grad[i] = _grad[cons_id_i];
    for (unsigned int j = 0; j < params.optimization.constraints.size(); ++j) {
      unsigned int cons_id_j = params.optimization.constraints[j].id;
      hessian[i][j] = _hess[cons_id_i][cons_id_j];
    }
  }

  // Factorise
  const unsigned int m = params.optimization.constraints.size();
  for (unsigned int ss = 0; ss < m - 1; ss++) {
    for (unsigned int i = ss + 1; i < m; i++) {
      hessian[i][ss] = hessian[i][ss] / hessian[ss][ss];
      for (unsigned int j = ss + 1; j < m; j++) {
        hessian[i][j] = hessian[i][j] - hessian[i][ss] * hessian[ss][j];
      }
    }
  }

  // Solve
  for (unsigned int i = 1; i < m; i++) {
    double a = 0.0;
    for (unsigned int j = 0; j < i; j++) {
      a = a - hessian[i][j] * grad[j];
    }
    grad[i] = grad[i] + a;
  }
  grad[m - 1] = grad[m - 1] / hessian[(m - 1)][(m - 1)];
  for (int i = m - 2; i >= 0; i--) {
    double a = grad[i];
    for (int j = i + 1; j < (int)m; j++) {
      a = a - hessian[i][j] * grad[j];
    }
    grad[i] = a / hessian[i][i];
  }

  unsigned int ii = 0;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id_i = params.optimization.constraints[i].id;
    _grad[cons_id_i] = grad[ii];
    ++ii;
  }

} // end of _Solve()

void MMA::_DualLineSearch() {
  auto params = parameters();
  double theta = 1.005;
  const unsigned int m = params.optimization.constraints.size() *
                         1000; // this varible has no significance this is just for indexing
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    if (theta < -1.01 * _s[cons_id] / _lam[cons_id]) {
      theta = -1.01 * _s[cons_id] / _lam[cons_id];
    }
    if (theta < -1.01 * _s[cons_id + m] / _mu[cons_id]) {
      theta = -1.01 * _s[cons_id + m] / _mu[cons_id];
    }
  }
  theta = 1.0 / theta;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    _lam[cons_id] = _lam[cons_id] + theta * _s[cons_id];
    _mu[cons_id] = _mu[cons_id] + theta * _s[cons_id + m];
  }
} // end of _DualLineSearch()

double MMA::_DualResidual() {
  auto params = parameters();
  const unsigned int m = params.optimization.constraints.size() *
                         1000; // this varible has no significance this is just for indexing
  std::map<unsigned int, double> res;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    res[cons_id] = 0.0;
    res[cons_id + m] = 0.0;
    for (const auto &[cell_id, value] : _pseudo_densities) {
      res[cons_id] += _pij[cons_id][cell_id] / (_upp_asym[cell_id] - _pseudo_densities[cell_id]) +
                      _qij[cons_id][cell_id] / (_pseudo_densities[cell_id] - _low_asym[cell_id]);
    }
  }

  // all reduce
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    res[cons_id] = Utilities::MPI::sum(res[cons_id], mpi_communicator());
  }

  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    res[cons_id] += -_b[cons_id] - _a[cons_id] * _z - _y[cons_id] + _mu[cons_id];
    res[cons_id + m] += _mu[cons_id] * _lam[cons_id] - _epsi;
  }
  double nrI = 0.0;
  for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
    unsigned int cons_id = params.optimization.constraints[i].id;
    if (nrI < std::abs(res[cons_id])) {
      nrI = std::abs(res[cons_id]);
    }
    if (nrI < std::abs(res[cons_id + m])) {
      nrI = std::abs(res[cons_id + m]);
    }
  }
  return nrI;
}