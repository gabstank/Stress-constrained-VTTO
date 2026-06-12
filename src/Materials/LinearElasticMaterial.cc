// C++ headers

// Deal.II headers
#include <deal.II/physics/elasticity/standard_tensors.h>

// Project headers
#include <DesignField.h>
#include <LinearElasticMaterial.h>
#include <MacrosAndTypedefs.h>
#include <Parameter.h>

MaterialLinElastic::MaterialLinElastic()
    : _lambda(parameters().material.elastic_isotropic.lambda),
      _mu(parameters().material.elastic_isotropic.mu) {

  _elasticity_tensor = _lambda * StandardTensors::IxI() + 2 * _mu * StandardTensors::II();
}

SymmetricTensor<2, 2> MaterialLinElastic::GetCauchyStress(const SymmetricTensor<2, 2> &e,
                                                          const double &pseudo_density) {
  return DesignField::GetPhysicalDensity(pseudo_density) * _elasticity_tensor * e;
}

double MaterialLinElastic::GetVonMisesStress(const SymmetricTensor<2, 2> &e,
                                             const double &pseudo_density,
                                             const double &stress_allow) {
  SymmetricTensor<2, 2> stress = GetCauchyStress(e, pseudo_density);
  // Stress_allow is used to ensure numerical stability of sensitivity analysis
  double vm_stress = std::sqrt(1.5 * stress * stress - 0.5 * std::pow(trace(stress), 2) +
                               1.e-4 * stress_allow * stress_allow);
  if (vm_stress < 0.0)
    throw std::runtime_error("von Mises stress should be always positive. Check in "
                             "MaterialLinElastic::GetVonMisesStress");
  return vm_stress;
}

double MaterialLinElastic::GetVonMisesStressFromCauchy(const SymmetricTensor<2, 2> &cauchy,
                                                       const double &stress_allow) {
  // Stress_allow is used to ensure numerical stability of sensitivity analysis
  double vm_stress = std::sqrt(1.5 * cauchy * cauchy - 0.5 * std::pow(trace(cauchy), 2) +
                               1.e-4 * stress_allow * stress_allow);
  if (vm_stress < 0.0)
    throw std::runtime_error("von Mises stress should be always positive. Check in "
                             "MaterialLinElastic::GetVonMisesStress");
  return vm_stress;
}