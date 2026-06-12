// Deal.II headers
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

// Project headers
#include <Filter.h>

class PDEFilter : public Filter {
public:
  PDEFilter(Mesh &, std::string);
  virtual ~PDEFilter() = default;
  virtual void updateContinuationParameters() override;
  virtual void applyFilter(std::map<CellId, double> &pseudo_densities) override;
  virtual void applyChainrule(std::map<CellId, double> &sensitivities) override;
  std::string getFilterName() const override {
    return FilterNames::PDEFilter;
  } // Override to return name
protected:
  virtual void _SetActiveIndex();
  void _SetupFilterSystem();
  void _AssembleFilterSystem();
  void _AssembleFilterSystemBoundary();
  void _AssembleFilterRHS();
  void _SolveFilterSystem();
  void _InterpolateCellValuesFromNodes(bool set_negative_values_to_zero = false);

  std::string _type;
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
  double _filter_decrease_rate;
};