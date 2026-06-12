#pragma once

// Deal.II headers
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/affine_constraints.h>

// Project headers
#include <DataOutput.h>
#include <Mesh.h>
#include <Miscellaneous.h>

class DataTransfer {
public:
  DataTransfer(Mesh &mesh_);
  virtual ~DataTransfer() = default;
  // Create a new DataUnit if necessary, initialize the sol_transfer
  void RegisterDataMap(const std::map<CellId, double> &data_map_in, const std::string &key);
  // check if key exists and has been transfered
  std::map<CellId, double> GetTransferedMap(const std::string &key);
  void PrepareForTransfer();
  // Reinit the dof handler, locally owned dofs etc...; interpolate; store in
  // the maps
  void ExecuteTransfer();
  void GetOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output);

private:
  DoFHandler<2> _density_dof_handler;
  IndexSet _locally_owned_dofs, _locally_relevant_dofs;
  FESystem<2> _fe_density;
  AffineConstraints<double> _density_constraints;

  struct DataUnit {
    DataUnit(const std::map<CellId, double> &in_map_, const DoFHandler<2> &density_dof_handler_)
        : in_out_map(in_map_), transfered(false), sol_transfer(density_dof_handler_) {}
    std::map<CellId, double> in_out_map;
    bool transfered;
    TrilinosWrappers::MPI::Vector in_out_vec;
    TrilinosWrappers::MPI::Vector transfer_vec;
    TrilinosWrappers::MPI::Vector vec_for_grad;
    parallel::distributed::SolutionTransfer<2, TrilinosWrappers::MPI::Vector> sol_transfer;
  };

  std::map<std::string, std::unique_ptr<DataUnit>> _transfer_data;
};