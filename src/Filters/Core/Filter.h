#ifndef FILTER_H
#define FILTER_H

// Project headers
#include <Mesh.h>
#include <Miscellaneous.h>
#include <Parameter.h>

class Filter {
public:
  Filter(Mesh &mesh_) : _mesh(mesh_) {}
  virtual ~Filter() = default;
  virtual void updateContinuationParameters() = 0;
  virtual void applyFilter(std::map<CellId, double> &pseudo_densities) = 0;
  virtual void applyChainrule(std::map<CellId, double> &sensitivities) = 0;
  virtual std::string getFilterName() const = 0; // Pure virtual function to get filter name
  void SetNonDesignCells(std::vector<CellId> non_design_cells) {
    _non_design_cells = non_design_cells;
  }

protected:
  Mesh &_mesh;
  std::vector<CellId> _non_design_cells;
};

#endif