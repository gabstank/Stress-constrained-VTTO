#pragma once

// C++ headers
#include <iostream>

// Deal.II headers
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/numerics/data_postprocessor.h>

#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/transformations.h>

// Project headers
#include <DataOutput.h>
#include <Mesh.h>
#include <Parameter.h>
#include <StressFilter.h>

class BVP {
public:
  BVP(Mesh &);

  virtual ~BVP() = default;

  virtual void Run() = 0;

  //-------------------------------------------------------------

  /**
   * @brief Function to get the base output data that are displayed in
   * optimization.
   *
   * @param output
   */
  virtual void GetAllOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output);

  /// TODO: make this Pure virtual function.
  virtual Vector<double> GetPostprocessingData(unsigned int &data_flag) {
    (void)data_flag;
    throw std::runtime_error("BVP::GetPostprocessingData should never be called in base class. "
                             "TODO: Once this funcion is implemented in all the classes this pure "
                             "virtual.");
  }

  virtual Vector<double> GetRefinementData(const std::string &ref_type);

public:
  void _SolveLinearSystem(LA::MPI::SparseMatrix &system_matrix, Vector<double> &solution,
                          LA::MPI::Vector &system_rhs);
  void _SetupSmoothingSystem();
  void AssembleSmoothingMatrix();
  void SolveSmoothingSystem();

  IndexSet _locally_owned_dofs, _locally_relevant_dofs;
  hp::FECollection<2> _fe_collection;
  hp::QCollection<2> _q_collection;
  hp::QCollection<1> _face_q_collection;

  //-----------------Postprocessing data-----------------
  // System of equations to smooth the strain and stress output
  IndexSet _locally_owned_smooth_dofs, _locally_relevant_smooth_dofs;
  FESystem<2> _smoothing_fe;
  hp::FECollection<2> _smoothing_fe_collection;
  hp::QCollection<2> _smoothing_q_collection;
  AffineConstraints<double> _smoothing_constraints;
  LA::MPI::SparseMatrix _smoothing_matrix;
  LA::MPI::Vector _smoothing_solution;
  LA::MPI::Vector _smoothing_rhs;
  DoFHandler<2> _smoothing_dof_handler;
  //-----------------Postprocessing data-----------------

  DomainParallelTriaType &_tria;
  DoFHandler<2> _dof_handler;

  //-----------------Postprocessing data for stresses-----------------
  virtual void
  _SetElementStresses() = 0; // Pure virtual function to be implemented in derived classes
  StressFilter _stress_filter;
  std::map<CellId, double> _el_stress_map_raw;
  std::map<CellId, double> _el_stress_map_smoothed;
  virtual std::map<CellId, double> GetElementStresses() { return _el_stress_map_smoothed; }
  //---------------------------------------------------

  AffineConstraints<double> _constraints;
  Vector<double> _solution;

  LA::MPI::SparseMatrix _system_matrix;
  LA::MPI::Vector _system_rhs;
  LA::MPI::Vector _external_force;

}; // end of BVP class definition