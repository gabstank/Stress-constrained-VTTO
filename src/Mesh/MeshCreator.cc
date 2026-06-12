#include <Mesh.h>
#include <MeshCreator.h>

void MeshCreator::_CreateHyperRectangleDomain(DomainParallelTriaTypeDistributed &tria) {

  auto params = parameters();
  if (params.geometry.padding) {
    _CreatePaddedHyperRectangleDomain(tria);
    return;
  }

  // Gather data from .prm file to Generate grid
  Point<2> p1;
  Point<2> p2;
  std::vector<unsigned int> subdivisions;

  p1 = Point<2>(params.geometry.hyper_rectangle.point1[0],
                params.geometry.hyper_rectangle.point1[1]);
  p2 = Point<2>(params.geometry.hyper_rectangle.point2[0],
                params.geometry.hyper_rectangle.point2[1]);
  subdivisions.push_back(params.geometry.hyper_rectangle.subdivision[0]);
  subdivisions.push_back(params.geometry.hyper_rectangle.subdivision[1]);

  tria.clear(); // Just to be sure I clear here, if deemed not useful can be
                // removed.

  GridGenerator::subdivided_hyper_rectangle(tria, subdivisions, p1, p2,
                                            false); // Generate the hyper_rectangle mesh.
}

void MeshCreator::_CreateHyperLDomain(DomainParallelTriaTypeDistributed &tria) {
  auto params = parameters();
  if (params.geometry.padding) {
    _CreatePaddedHyperLDomain(tria);
    return;
  }

  // Gather data from .prm file to Generate grid
  Point<2> p1;
  Point<2> p2;
  std::vector<unsigned int> subdivisions;
  std::vector<int> cells_to_remove;

  p1 = Point<2>(params.geometry.hyper_l.point1[0], params.geometry.hyper_l.point1[1]);
  p2 = Point<2>(params.geometry.hyper_l.point2[0], params.geometry.hyper_l.point2[1]);
  subdivisions.push_back(params.geometry.hyper_l.subdivision[0]);
  subdivisions.push_back(params.geometry.hyper_l.subdivision[1]);

  cells_to_remove.push_back(params.geometry.hyper_l.cells_to_remove[0]);
  cells_to_remove.push_back(params.geometry.hyper_l.cells_to_remove[1]);

  tria.clear(); // Just to be sure I clear here, if deemed not useful can be
                // removed.

  GridGenerator::subdivided_hyper_L(tria, subdivisions, p1, p2,
                                    cells_to_remove); // Generate the hyper_L mesh.
}

void MeshCreator::_CreateAbaqusDomain(DomainParallelTriaTypeDistributed &tria) {

  auto params = parameters();

  GridIn<2> p_grid_in;
  std::ifstream p_input_file(
      params.geometry.abaqus.input); // rereading the file because of dealii I/O error.
                                     // Could not use the input_file again.
  p_grid_in.attach_triangulation(tria);
  p_grid_in.read_abaqus(p_input_file);
}

void MeshCreator::_RemovePaddedCellsForBoundaryConditions(DomainParallelTriaTypeDistributed &tria) {
  auto params = parameters();

  std::map<int, int> bid_counter;
  for (unsigned int i = 0; i < params.geometry.boundary_ids.size(); ++i) {
    if (bid_counter.find(params.geometry.boundary_ids[i]) == bid_counter.end())
      bid_counter.insert(std::pair<int, int>(params.geometry.boundary_ids[i], 2));
    else
      bid_counter.find(params.geometry.boundary_ids[i])->second -= 1;
  }
  std::set<typename DomainParallelTriaTypeDistributed::active_cell_iterator> cells_to_remove_set;
  for (auto cell : tria.active_cell_iterators()) {
    if (cell->at_boundary())
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        int boundary_id = Mesh::GetBoundaryID(bid_counter, cell->face(f)->center());
        if (boundary_id > 0) {
          cells_to_remove_set.insert(cell);
          break;
        }
      }
  }

  GridGenerator::create_triangulation_with_removed_cells(tria, cells_to_remove_set,
                                                         tria); // Generate the hyper_L mesh.
}

void MeshCreator::_CreatePaddedHyperRectangleDomain(DomainParallelTriaTypeDistributed &tria) {
  auto params = parameters();

  // Gather data from .prm file to Generate grid
  Point<2> p1;
  Point<2> p2;
  std::vector<unsigned int> subdivisions;

  double cell_size =
      (params.geometry.hyper_rectangle.point2[0] - params.geometry.hyper_rectangle.point1[0]) /
      params.geometry.hyper_rectangle.subdivision[0];
  p1 = Point<2>(params.geometry.hyper_rectangle.point1[0] - cell_size,
                params.geometry.hyper_rectangle.point1[1] - cell_size);
  p2 = Point<2>(params.geometry.hyper_rectangle.point2[0] + cell_size,
                params.geometry.hyper_rectangle.point2[1] + cell_size);
  subdivisions.push_back(params.geometry.hyper_rectangle.subdivision[0] + 2);
  subdivisions.push_back(params.geometry.hyper_rectangle.subdivision[1] + 2);

  tria.clear(); // Just to be sure I clear here, if deemed not useful can be
                // removed.

  GridGenerator::subdivided_hyper_rectangle(tria, subdivisions, p1, p2,
                                            false); // Generate the hyper_rectangle mesh.

  // Assign material ID for the cells in the padding region to be able to identify them later
  for (auto cell : tria.active_cell_iterators()) {
    if (cell->center()[0] < params.geometry.hyper_rectangle.point1[0] ||
        cell->center()[0] > params.geometry.hyper_rectangle.point2[0] ||
        cell->center()[1] < params.geometry.hyper_rectangle.point1[1] ||
        cell->center()[1] > params.geometry.hyper_rectangle.point2[1]) {
      cell->set_material_id(e_padding_mat_id); // Marking the cells outside the original domain with
                                               // a different material ID for now.
    } else
      cell->set_material_id(0);
  }

  // Remove the cells for the cutout to avoid issues with the filter and the boundary conditions.
  _RemovePaddedCellsForBoundaryConditions(tria);
}

void MeshCreator::_CreatePaddedHyperLDomain(DomainParallelTriaTypeDistributed &tria) {
  auto params = parameters();

  // Gather data from .prm file to Generate grid
  Point<2> p1;
  Point<2> p2;
  std::vector<unsigned int> subdivisions;
  std::vector<int> cells_to_remove;

  double cell_size = (params.geometry.hyper_l.point2[0] - params.geometry.hyper_l.point1[0]) /
                     params.geometry.hyper_l.subdivision[0];
  p1 = Point<2>(params.geometry.hyper_l.point1[0] - cell_size,
                params.geometry.hyper_l.point1[1] - cell_size);
  p2 = Point<2>(params.geometry.hyper_l.point2[0] + cell_size,
                params.geometry.hyper_l.point2[1] + cell_size);
  subdivisions.push_back(params.geometry.hyper_l.subdivision[0] + 2);
  subdivisions.push_back(params.geometry.hyper_l.subdivision[1] + 2);
  cells_to_remove.push_back(params.geometry.hyper_l.cells_to_remove[0]);
  cells_to_remove.push_back(params.geometry.hyper_l.cells_to_remove[1]);

  tria.clear(); // Just to be sure I clear here, if deemed not useful can be
                // removed.

  GridGenerator::subdivided_hyper_L(tria, subdivisions, p1, p2,
                                    cells_to_remove); // Generate the hyper_L mesh.

  // Assign material ID for the cells in the padding region to be able to identify them later
  double l_cutout_x = -params.geometry.hyper_l.cells_to_remove[0] *
                      ((params.geometry.hyper_l.point2[0] - params.geometry.hyper_l.point1[0]) /
                       params.geometry.hyper_l.subdivision[0]);
  double l_cutout_y = -params.geometry.hyper_l.cells_to_remove[1] *
                      ((params.geometry.hyper_l.point2[1] - params.geometry.hyper_l.point1[1]) /
                       params.geometry.hyper_l.subdivision[1]);

  for (auto cell : tria.active_cell_iterators()) {
    if (cell->center()[0] < params.geometry.hyper_l.point1[0] ||
        cell->center()[0] > params.geometry.hyper_l.point2[0] ||
        cell->center()[1] < params.geometry.hyper_l.point1[1] ||
        cell->center()[1] > params.geometry.hyper_l.point2[1]) {
      cell->set_material_id(e_padding_mat_id); // Marking the cells outside the original domain with
                                               // a different material ID for now.
    }
    // Consider the padding for the domain edges after removing the cells for the L cutout, to avoid
    // issues with the filter and the boundary conditions.

    else if (std::abs(cell->center()[0] - params.geometry.hyper_l.point2[0]) < l_cutout_x &&
             std::abs(cell->center()[1] - params.geometry.hyper_l.point2[1]) < l_cutout_y) {
      cell->set_material_id(e_padding_mat_id);
    } else
      cell->set_material_id(0);
  }

  // Remove the cells for the L cutout to avoid issues with the filter and the boundary conditions.
  _RemovePaddedCellsForBoundaryConditions(tria);
}
