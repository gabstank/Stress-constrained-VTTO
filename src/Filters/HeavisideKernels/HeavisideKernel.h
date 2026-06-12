#ifndef HEAVISIDEKERNEL_H
#define HEAVISIDEKERNEL_H

#include <map>
#include <vector>

// Deal.II headers
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

// Project Headers
#include <Mesh.h>

class HeavisideKernel {
public:
  virtual ~HeavisideKernel() = default;

  // Pure virtual functions for the specific math
  virtual void Apply(std::map<CellId, double> &density_map,
                     const std::map<CellId, double> &input_map, double beta,
                     const std::vector<CellId> &non_design_cells, Mesh &mesh) = 0;

  virtual void ChainRule(std::map<CellId, double> &sensitivity_map,
                         const std::map<CellId, double> &input_sensitivities,
                         const std::map<CellId, double> &physical_densities, double beta,
                         const std::vector<CellId> &non_design_cells, Mesh &mesh) = 0;
};

#endif