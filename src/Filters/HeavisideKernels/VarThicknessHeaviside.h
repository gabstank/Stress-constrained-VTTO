#include <HeavisideKernel.h>

class VarThicknessHeaviside : public HeavisideKernel {
public:
  virtual ~VarThicknessHeaviside() = default;

  // Pure virtual functions for the specific math
  virtual void Apply(std::map<CellId, double> &density_map,
                     const std::map<CellId, double> &input_map, double beta,
                     const std::vector<CellId> &non_design_cells, Mesh &mesh) override;

  virtual void ChainRule(std::map<CellId, double> &sensitivity_map,
                         const std::map<CellId, double> &input_sensitivities,
                         const std::map<CellId, double> &physical_densities, double beta,
                         const std::vector<CellId> &non_design_cells, Mesh &mesh) override;
};