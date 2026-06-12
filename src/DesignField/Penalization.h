#ifndef PENALIZATION_H
#define PENALIZATION_H

class Penalization {
public:
  virtual ~Penalization() = default;

  // Updates the global penalty parameter based on the strategy
  // We pass 'current_penalty' by reference so the strategy can modify it
  virtual void InitializePenalty(double &current_penalty) = 0;
  virtual void UpdatePenalty(double &current_penalty, const double design_change) = 0;
  virtual void InitializeStressExponent(double &current_exponent) = 0;
  virtual void UpdateStressExponent(double &current_exponent) = 0;

  virtual double GetPhysicalDensity(const double density, const double penalty) = 0;
  virtual double GetPhysicalDensityDer(const double density, const double penalty) = 0;
  virtual double GetStressRelaxation(const double density, const double penalty) = 0;
  virtual double GetStressRelaxationDer(const double density, const double penalty) = 0;

  virtual bool isSolidCell(const double density) = 0;
};

#endif // PENALIZATION_H