// C++ headers

// Deal.II headers
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/hp/refinement.h>

// Project headers
#include <DesignField.h>
#include <Mesh.h>
#include <MeshCreator.h>

// Implementation of static member variable (MUST be outside the class definition)
std::map<unsigned int, unsigned int> Mesh::symmetry_id_to_symmetry_plane;

Mesh::Mesh()
    : tria(mpi_communicator(),
           typename Triangulation<2>::MeshSmoothing(Triangulation<2>::smoothing_on_refinement |
                                                    Triangulation<2>::smoothing_on_coarsening),
           DomainParallelTriaTypeDistributed::default_setting
           /*DomainParallelTriaTypeDistributed::no_automatic_repartitioning*/) {

  // Generate mesh and output
  SetupMesh();
}

// destructor
Mesh::~Mesh() {}

void Mesh::SetupMesh() {

  MeshCreator::createMesh(tria, parameters().geometry.name);

  tria.refine_global(parameters().geometry.global_refinements);

  SetupBoundaryIDs();

  // printing mesh statistics
  if (parameters().debug.verbose) {
    std::vector<unsigned int> domain_cells_info =
        Utilities::MPI::gather(mpi_communicator(), tria.n_active_cells());
    std::vector<unsigned int> domain_vertices_info =
        Utilities::MPI::gather(mpi_communicator(), tria.n_vertices());

    pcout() << std::endl;
    pcout() << "Domain Parallel: No of cells: ";
    for (const auto &cell_info : domain_cells_info)
      pcout() << cell_info << " , ";
    pcout() << " | No of vertices: ";
    for (const auto &vert_info : domain_vertices_info)
      pcout() << vert_info << " , ";
    pcout() << std::endl;
    PROJ_MPI_BARRIER
  }
}

void Mesh::SetupBoundaryIDs() {
  // Mark the boundaries.
  // Defining some help variable
  // Map of {ID, counter}
  /**
   * The reason for having this map.
   * The counter is initial set to 2 for each unique ID.
   * If there are multiple conditions for an ID, each additional condition the
   * counter is reduced by 1.
   *
   * During the check if a condition is satisfied the counter will be increased
   * by 1 for corresponding ID. Finally whichever ID has counter value 2 + 1
   * is the Boundary ID of that point.
   */
  auto params = parameters();
  AssertThrow(params.geometry.boundary_ids.size() == params.geometry.boundary_comp.size(),
              ExcDimensionMismatch(params.geometry.boundary_ids.size(),
                                   params.geometry.boundary_comp.size()));
  AssertThrow(params.geometry.boundary_ids.size() == params.geometry.boundary_coord.size(),
              ExcDimensionMismatch(params.geometry.boundary_ids.size(),
                                   params.geometry.boundary_coord.size()));
  AssertThrow(params.geometry.boundary_ids.size() == params.geometry.boundary_tol.size(),
              ExcDimensionMismatch(params.geometry.boundary_ids.size(),
                                   params.geometry.boundary_tol.size()));

  std::map<int, int> bid_counter;
  for (unsigned int i = 0; i < params.geometry.boundary_ids.size(); ++i) {
    if (bid_counter.find(params.geometry.boundary_ids[i]) == bid_counter.end())
      bid_counter.insert(std::pair<int, int>(params.geometry.boundary_ids[i], 2));
    else
      bid_counter.find(params.geometry.boundary_ids[i])->second -= 1;
  }

  SetupSymmetryIDs();

  // For parallel domain
  // Loop over the cells. If the cell is at boundary, get boundary ID from the
  // function GetBoundaryID(...).
  Vector<double> p_dom_boundary_ids(tria.n_active_cells() * GeometryInfo<2>::faces_per_cell);
  Vector<double> p_out_dom_boundary_ids(tria.n_active_cells());
  unsigned int counter = 0, cell_counter = 0;
  for (auto cell : tria.active_cell_iterators()) {
    //      if(!cell->is_locally_owned()){
    //        ++cell_counter;
    //        counter += GeometryInfo<2>::faces_per_cell;
    //        continue;
    //      }
    if (cell->at_boundary() &&
        cell->material_id() != e_padding_mat_id) // Only mark the boundary cells that are not in the
                                                 // padding region
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        if (cell->face(f)->at_boundary()) {
          int boundary_id = GetBoundaryID(bid_counter, cell->face(f)->center());
          cell->face(f)->set_boundary_id(boundary_id);
          if (p_dom_boundary_ids[counter] == 0) {
            p_dom_boundary_ids[counter] = boundary_id;
            p_out_dom_boundary_ids[cell_counter] = boundary_id;
          }
        }
        ++counter;
      }
    ++cell_counter;
  }
}

int Mesh::GetBoundaryID(std::map<int, int> bid_counter, const Point<2> &point) {
  auto params = parameters();
  for (unsigned int i = 0; i < params.geometry.boundary_ids.size(); ++i) {
    if (std::abs(point[params.geometry.boundary_comp[i]] - params.geometry.boundary_coord[i]) <
        params.geometry.boundary_tol[i]) {
      bid_counter.find(params.geometry.boundary_ids[i])->second++;
      if (bid_counter.find(params.geometry.boundary_ids[i])->second == 2 + 1) {
        return bid_counter.find(params.geometry.boundary_ids[i])->first;
      }
    }
  }

  return 0;
}

void Mesh::SetupSymmetryIDs() {
  auto params = parameters();
  for (unsigned int i = 0; i < params.bvp.dirichlet.id.size(); ++i) {
    if (i == 0) {
      if (params.bvp.dirichlet.id[i] == params.bvp.dirichlet.id[i + 1])
        continue;
    } else if (i == params.bvp.dirichlet.id.size() - 1) {
      if (params.bvp.dirichlet.id[i] == params.bvp.dirichlet.id[i - 1])
        continue;
    } else {
      if (params.bvp.dirichlet.id[i] == params.bvp.dirichlet.id[i - 1] ||
          params.bvp.dirichlet.id[i] == params.bvp.dirichlet.id[i + 1])
        continue;
    }
    symmetry_id_to_symmetry_plane[params.bvp.dirichlet.id[i]] = params.bvp.dirichlet.component[i];
  }
}

void Mesh::ImposeSymmetryConditions(LA::MPI::Vector &vector, const DoFHandler<2> &dof_handler) {
  auto params = parameters();
  std::vector<bool> dof_touched(dof_handler.n_dofs());
  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    // if cell belong to some other domain, do nothing
    if (!cell->is_locally_owned())
      continue;
    if (cell->at_boundary())
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
        if (cell->face(f)->at_boundary() &&
            symmetry_id_to_symmetry_plane.find(cell->face(f)->boundary_id()) !=
                symmetry_id_to_symmetry_plane.end()) {
          std::vector<types::global_dof_index> face_dof_indices;
          face_dof_indices.resize(cell->get_fe().dofs_per_face);
          cell->face(f)->get_dof_indices(face_dof_indices, cell->active_fe_index());
          for (auto dof_index : face_dof_indices) {
            if (dof_touched[dof_index] || !locally_owned_dofs.is_element(dof_index))
              continue;
            vector[dof_index] *= 2;
            dof_touched[dof_index] = true;
          }
        }
      }
  }
}

void Mesh::hAdaptationCNF(const Vector<double> &vec_field_cnf_cnf,
                          const DoFHandler<2> &dof_handler) {
  auto params = parameters();
  auto densities = DesignField::getCurrent();

  std::vector<bool> vertices_for_ref(tria.n_vertices());
  std::vector<bool> vertices_for_coarsen(tria.n_vertices());
  double linfty_norm_local = 0, linfty_norm = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      Tensor<1, 2> config_force_v;
      for (unsigned int d = 0; d < 2; ++d) {
        unsigned int dof_index = cell->vertex_dof_index(v, d, cell->active_fe_index());
        config_force_v[d] = vec_field_cnf_cnf[dof_index];
      }
      if (config_force_v.norm() > linfty_norm_local)
        linfty_norm_local = config_force_v.norm();
    }
  }
  linfty_norm = Utilities::MPI::max(linfty_norm_local, mpi_communicator());

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      Tensor<1, 2> config_force_v;
      for (unsigned int d = 0; d < 2; ++d) {
        unsigned int dof_index = cell->vertex_dof_index(v, d, cell->active_fe_index());
        config_force_v[d] = vec_field_cnf_cnf[dof_index];
      }
      if (config_force_v.norm() >
              linfty_norm * params.optimization.mesh_adaptivity.relative_threshold_refine &&
          !vertices_for_ref[cell->vertex_index(v)]) {
        vertices_for_ref[cell->vertex_index(v)] = true;
      } else if (config_force_v.norm() <
                     linfty_norm * params.optimization.mesh_adaptivity.relative_threshold_coarsen &&
                 (densities.at(cell->id()) < 0.01 || densities.at(cell->id()) > 0.99) &&
                 !vertices_for_coarsen[cell->vertex_index(v)]) {
        vertices_for_coarsen[cell->vertex_index(v)] = true;
      }
    }
  }
  for (auto cell : tria.active_cell_iterators()) {
    if (cell->is_artificial())
      continue;
    bool all_vertices_for_coarsen = true;
    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      if (vertices_for_ref[cell->vertex_index(v)] && !cell->refine_flag_set() &&
          cell->level() < (int)(params.optimization.mesh_adaptivity.max_h_refinements))
        cell->set_refine_flag();
      if (!vertices_for_coarsen[cell->vertex_index(v)])
        all_vertices_for_coarsen = false;
    }
    if (all_vertices_for_coarsen && !cell->coarsen_flag_set() && !cell->refine_flag_set())
      cell->set_coarsen_flag();
  }
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  SetupBoundaryIDs();
}

void Mesh::hAdaptationVNM(const Vector<double> &scalar_field_vnm,
                          const DoFHandler<2> &scalar_dof_handler) {
  auto params = parameters();

  auto densities = DesignField::get(DesignFieldNames::FilteredDensity);

  std::vector<bool> vertices_for_ref(tria.n_vertices());
  std::vector<bool> vertices_for_coarsen(tria.n_vertices());
  double linfty_norm_local = 0, linfty_norm = 0;
  for (const auto &cell : scalar_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    double cell_density = densities.at(cell->id());
    double cell_stress_relaxation = DesignField::StressRelaxation(cell_density);
    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      unsigned int dof_index = cell->vertex_dof_index(v, 0, cell->active_fe_index());
      if (scalar_field_vnm[dof_index] * cell_stress_relaxation > linfty_norm_local)
        linfty_norm_local = cell_stress_relaxation * scalar_field_vnm[dof_index];
    }
  }
  linfty_norm = Utilities::MPI::max(linfty_norm_local, mpi_communicator());
  for (const auto &cell : scalar_dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    // double cell_density = densities.at(cell->id());
    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      unsigned int dof_index = cell->vertex_dof_index(v, 0, cell->active_fe_index());
      if (scalar_field_vnm[dof_index] >
              linfty_norm * params.optimization.mesh_adaptivity.relative_threshold_refine &&
          !vertices_for_ref[cell->vertex_index(v)]) {
        vertices_for_ref[cell->vertex_index(v)] = true;
      } else if (scalar_field_vnm[dof_index] <
                     linfty_norm * params.optimization.mesh_adaptivity.relative_threshold_coarsen &&
                 !vertices_for_coarsen[cell->vertex_index(v)]) {
        vertices_for_coarsen[cell->vertex_index(v)] = true;
      }
    }
  }
  for (auto cell : tria.active_cell_iterators()) {
    if (cell->is_artificial())
      continue;
    bool all_vertices_for_coarsen = true;
    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
      if (vertices_for_ref[cell->vertex_index(v)] && !cell->refine_flag_set() &&
          cell->level() < (int)(params.optimization.mesh_adaptivity.max_h_refinements))
        cell->set_refine_flag();
      if (!vertices_for_coarsen[cell->vertex_index(v)])
        all_vertices_for_coarsen = false;
    }
    if (all_vertices_for_coarsen && !cell->coarsen_flag_set() && !cell->refine_flag_set())
      cell->set_coarsen_flag();
  }
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  SetupBoundaryIDs();
}

void Mesh::hAdaptationDENS() {
  auto params = parameters();

  auto densities = DesignField::get(DesignFieldNames::FilteredDensity);

  for (auto cell : tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;
    if (densities.at(cell->id()) <= 0.8 && densities.at(cell->id()) >= 0.2 &&
        cell->level() < (int)(params.optimization.mesh_adaptivity.max_h_refinements))
      cell->set_refine_flag();
    else if (densities.at(cell->id()) <= 0.01 || densities.at(cell->id()) >= 0.99)
      cell->set_coarsen_flag();
  }
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  SetupBoundaryIDs();
}

void Mesh::hAdaptationDensityJump() {
  auto params = parameters();

  auto densities_local = DesignField::getCurrent();

  // reduce density data
  std::vector<std::map<CellId, double>> temp_input =
      Utilities::MPI::all_gather(mpi_communicator(), densities_local);
  std::map<CellId, double> densities;
  for (const auto &input_map : temp_input)
    for (const auto &[cellid, val] : input_map)
      if (densities.find(cellid) == densities.end())
        densities[cellid] = val;

  Vector<double> max_density_jump_per_cell(tria.n_active_cells());

  for (auto cell : tria.active_cell_iterators()) {
    if (!cell->is_locally_owned())
      continue;

    for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
      if (cell->face(f)->at_boundary())
        continue;

      auto neighbor = cell->neighbor(f);

      if (neighbor->has_children()) { // check if neighbor has higher refinement level
        for (unsigned int neighbor_child_id = 0; neighbor_child_id < neighbor->n_children();
             ++neighbor_child_id) {
          for (unsigned int f_neighbor = 0; f_neighbor < GeometryInfo<2>::faces_per_cell;
               ++f_neighbor) {
            if (neighbor->child(neighbor_child_id)->face(f_neighbor)->at_boundary())
              continue;
            if (neighbor->child(neighbor_child_id)->neighbor(f_neighbor)->id() ==
                cell->id()) { // find the children of the neighbor that are actually neighboring the
                              // current cell
              double jump = std::abs(densities[cell->id()] -
                                     densities[neighbor->child(neighbor_child_id)->id()]);
              if (max_density_jump_per_cell[cell->active_cell_index()] < jump)
                max_density_jump_per_cell[cell->active_cell_index()] = jump;
            }
          }
        }
      } else {
        double jump = std::abs(densities[cell->id()] - densities[neighbor->id()]);
        if (max_density_jump_per_cell[cell->active_cell_index()] < jump)
          max_density_jump_per_cell[cell->active_cell_index()] = jump;
      }
    }
  }

  unsigned int n_levels = 1;
  GridRefinement::refine(tria, max_density_jump_per_cell,
                         params.optimization.mesh_adaptivity.relative_threshold_refine /
                             n_levels); // Refine threshold adjusted to the number of thicknesses
  GridRefinement::coarsen(tria, max_density_jump_per_cell,
                          params.optimization.mesh_adaptivity.relative_threshold_coarsen /
                              n_levels); // Coarsen threshold adjusted to the number of thicknessess

  // Limit maximum h refinement
  Assert(tria.n_levels() <= params.optimization.mesh_adaptivity.max_h_refinements + 1,
         ExcInternalError());

  if (tria.n_levels() > params.optimization.mesh_adaptivity.max_h_refinements)
    for (const auto &cell :
         tria.active_cell_iterators_on_level(params.optimization.mesh_adaptivity.max_h_refinements))
      cell->clear_refine_flag(); // Limit h refinements from top
  for (const auto &cell :
       tria.active_cell_iterators_on_level(params.optimization.mesh_adaptivity.min_h_refinements))
    cell->clear_coarsen_flag(); // Limit h refinements from bottom

  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  SetupBoundaryIDs();
}

void Mesh::RefineAll() {
  auto params = parameters();

  tria.set_all_refine_flags();
  // Limit maximum h refinement
  Assert(tria.n_levels() <= params.optimization.mesh_adaptivity.max_h_refinements + 1,
         ExcInternalError());

  if (tria.n_levels() > params.optimization.mesh_adaptivity.max_h_refinements)
    for (const auto &cell :
         tria.active_cell_iterators_on_level(params.optimization.mesh_adaptivity.max_h_refinements))
      cell->clear_refine_flag(); // Limit h refinements

  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  PROJ_MPI_BARRIER
}