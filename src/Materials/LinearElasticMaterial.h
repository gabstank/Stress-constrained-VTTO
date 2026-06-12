#pragma once

// C++ headers

// Project headers
#include <Miscellaneous.h>

class MaterialLinElastic {
public:
  MaterialLinElastic();
  SymmetricTensor<4, 2> GetElasticityTensor() { return _elasticity_tensor; }
  SymmetricTensor<2, 2> GetCauchyStress(const SymmetricTensor<2, 2> &e,
                                        const double &pseudo_density = 1.0);
  double GetVonMisesStress(const SymmetricTensor<2, 2> &e, const double &pseudo_density = 1.0,
                           const double &stress_allow = 0.0);
  double GetVonMisesStressFromCauchy(const SymmetricTensor<2, 2> &cauchy,
                                     const double &vm_allow = 0.);

private:
  const double _lambda, _mu;
  SymmetricTensor<4, 2> _elasticity_tensor;
};