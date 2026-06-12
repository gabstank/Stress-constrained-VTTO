#include <Penalization.h>

class VariableThickness : public Penalization {
public:
  void InitializePenalty(double &current_penalty) override;
  void UpdatePenalty(double &current_penalty, const double) override;
  void InitializeStressExponent(double &current_exponent) override;
  void UpdateStressExponent(double &current_exponent) override;
  double GetPhysicalDensity(const double density, const double penalty) override;
  double GetPhysicalDensityDer(const double density, const double penalty) override;
  double GetStressRelaxation(const double density, const double penalty) override;
  double GetStressRelaxationDer(const double density, const double penalty) override;
  bool isSolidCell(const double density) override;

private:
  double _penalty_growth_rate = 1.; // Set in InitializePenalty if penalty continuation is active
  double _target_stress_relaxation_exponent;           // Set in InitializeStressExponent
  double _stress_relaxation_exponent_growth_rate = 1.; // Set in Initialize
};