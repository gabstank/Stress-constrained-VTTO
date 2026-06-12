#include <VarThicknessHeaviside.h>

void VarThicknessHeaviside::Apply(std::map<CellId, double> &density_map,
                                  const std::map<CellId, double> &input_map, double beta,
                                  const std::vector<CellId> &non_design_cells, Mesh &mesh) {
  if (beta < 1.)
    beta = 1.;
  double low_thickness =
      parameters().optimization.design_field.variable_thickness.low_thickness_threshold;
  for (const auto &cell : mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (std::find(non_design_cells.begin(), non_design_cells.end(), cell->id()) !=
        non_design_cells.end()) {
      density_map[cell->id()] = input_map.at(cell->id());
      continue;
    }

    double switch_function =
        0.5 * (1. + std::tanh(beta * (input_map.at(cell->id()) / low_thickness -
                                      std::pow(low_thickness, 1. / beta))));

    density_map[cell->id()] = (1. - switch_function) * std::pow(input_map.at(cell->id()), beta) +
                              switch_function * input_map.at(cell->id());
  } // end of loop over cells
}

void VarThicknessHeaviside::ChainRule(std::map<CellId, double> &sensitivity_map,
                                      const std::map<CellId, double> &input_sensitivities,
                                      const std::map<CellId, double> &physical_densities,
                                      double beta, const std::vector<CellId> &non_design_cells,
                                      Mesh &mesh) {
  auto params = parameters();

  double low_thickness =
      params.optimization.design_field.variable_thickness.low_thickness_threshold;
  for (const auto &cell : mesh.tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (std::find(non_design_cells.begin(), non_design_cells.end(), cell->id()) !=
        non_design_cells.end()) {
      sensitivity_map[cell->id()] = input_sensitivities.at(cell->id());
      continue;
    }

    double switch_function =
        0.5 * (1. + std::tanh(beta * (physical_densities.at(cell->id()) / low_thickness -
                                      std::pow(low_thickness, 1. / beta))));
    double switch_function_der =
        beta / (2. * low_thickness) *
        (1. - std::pow(std::tanh(beta * (physical_densities.at(cell->id()) / low_thickness -
                                         std::pow(low_thickness, 1. / beta))),
                       2));

    double projection_der =
        switch_function_der * (physical_densities.at(cell->id()) -
                               std::pow(physical_densities.at(cell->id()), beta)) +
        (1. - switch_function) * beta * std::pow(physical_densities.at(cell->id()), beta - 1.) +
        switch_function;

    sensitivity_map[cell->id()] = input_sensitivities.at(cell->id()) * projection_der;
  }
}
