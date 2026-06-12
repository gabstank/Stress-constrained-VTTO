#pragma once

// Project headers
#include <DataOutput.h>
#include <MacrosAndTypedefs.h>

class Mesh {
public:
  // Static methods that don't depend on instance state
  static int GetBoundaryID(std::map<int, int> bid_counter, const Point<2> &point);
  static void ImposeSymmetryConditions(LA::MPI::Vector &vector, const DoFHandler<2> &dof_handler);

  Mesh();
  ~Mesh();

  // Instance methods that operate on the mesh state
  // h adaptivity functionality
  void hAdaptationCNF(const Vector<double> &vec_field_cnf, const DoFHandler<2> &dof_handler);
  void hAdaptationVNM(const Vector<double> &scalar_field_vnm,
                      const DoFHandler<2> &scalar_dof_handler);
  void hAdaptationDENS();
  void hAdaptationDensityJump();

  void RefineAll();

  // Public member variable
  DomainParallelTriaTypeDistributed tria; // Parallel domain mesh

private:
  void SetupMesh();
  void SetupBoundaryIDs();

  // Static method and variable
  static void SetupSymmetryIDs();
  static std::map<unsigned int, unsigned int> symmetry_id_to_symmetry_plane;
};