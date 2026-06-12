// Deal.II headers
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

// Project headers
#include <Mesh.h>

class StressFilter {
public:
  StressFilter(DomainParallelTriaType &);
  ~StressFilter() = default;
  void applyFilter(std::map<CellId, double> &pseudo_densities);
  std::string getFilterName() const { return FilterNames::StressFilter; } // Override to return name
protected:
  void _SetActiveIndex();
  void _SetupFilterSystem();
  void _AssembleFilterSystem();
  void _AssembleFilterRHS();
  void _SolveFilterSystem();
  void _InterpolateCellValuesFromNodes();

  IndexSet _locally_owned_dofs, _locally_relevant_dofs;
  DoFHandler<2> _filter_dof_handler;
  // New
  hp::FECollection<2> _fe_collection;
  FESystem<2> _fe_nothing;
  // Standard
  FESystem<2> _fe_filter;
  AffineConstraints<double> _filter_constraints;
  SparsityPattern _sparsity_pattern;

  std::map<CellId, double> _in_out_map;
  LA::MPI::Vector _filtered_field;

  LA::MPI::SparseMatrix _filter_matrix;
  LA::MPI::Vector _filter_rhs;

  double _filter_length_parameter;
  bool _solid_cell_exists = false;
};