#ifndef DESIGNFIELD_H
#define DESIGNFIELD_H

#include <Miscellaneous.h>
#include <Penalization.h>

class DesignField {
public:
  // No constructor or destructor for a fully static class usually
  // Or make them private if they have a specific static initialization/cleanup role
  // Initialization method - MUST be called at start of program
  static void Initialize();

  static void UpdatePenalty(const double design_change);
  static double GetPhysicalDensity(const double density);
  static double GetPhysicalDensityDer(const double density);
  static double StressRelaxation(const double density);
  static double StressRelaxationDer(const double density);

  static bool isSolidCell(const double density);

  static std::map<CellId, double> get(const std::string &name);
  static std::map<CellId, double> getCurrent();
  static std::map<std::string, std::map<CellId, double>> getAllActive();
  inline static double getPenalty() { return _penalty; }
  static void set(const std::map<CellId, double> &design_field, const std::string &name);

private:
  // State
  static std::unique_ptr<Penalization> _model; // The Strategy
  static double _penalty;
  static double _stress_penalty;

  // Static member variables
  static std::map<CellId, double> _pseudo_densities;
  static std::map<CellId, double> _filtered_densities;
  static std::map<CellId, double> _vt_projected_densities;
  static std::map<CellId, double> _heaviside_projected_densities;
  static std::map<CellId, double> _prev_pseudo_densities;
};

#endif // DESIGNFIELD_H