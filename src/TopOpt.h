#pragma once

// C++ headers
#include <iostream>
#include <memory> // smart pointers

// Deal.II headers
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/grid/grid_refinement.h>

// Project headers
#include <BVP.h>
#include <BVPFactory.h>
#include <DataOutput.h>
#include <DataTransfer.h>
#include <DesignField.h>
#include <ElasticityLin.h>
#include <Filter.h>
#include <OptimizerFactory.h>
#include <Parameter.h>
#include <ResponseHandler.h>

class TopOpt {
public:
  TopOpt();
  ~TopOpt() = default;
  void Run();

private:
  //----------- Optimization preparation ------------------
  void _InitializeBVP(const std::string &bvp_name);
  void _SetNonDesignCells();
  void _SetInitialDensity();
  void _SetFilters();
  void _PrepareOptimization();
  //----------- Optimization iteration steps --------------
  void _ComputeResponseValuesGradients();
  void _applyChainrule();
  void _DesignUpdate();
  void _FilterDensities();
  void _ComputeStoppingCriterion();
  void _AdaptiveRefinement();
  void _UpdateContinuationParameters();
  // ---------- Output and logging functionality ----------
  void _PushMapdataToOutput(const std::map<CellId, double> &, std::string,
                            DataOutput<DomainParallelTriaType, DoFHandler<2>> &);
  void _WriteDomainOutput(const unsigned int &iter);
  void _AddDensityDataToOutput();
  void _WriteIterationData();
  // ----------- Post processing --------------------------
  void _PrepareForPostprocessing();

  Mesh _mesh;
  std::unique_ptr<BVP> _bvp;
  std::unique_ptr<Optimizer> _optimizer;

  std::vector<std::unique_ptr<Filter>> _filters;

  DataTransfer _data_transfer;
  ResponseHandler _response_handler;
  std::unique_ptr<DataOutput<DomainParallelTriaType, DoFHandler<2>>> _output;

  std::vector<CellId> _non_design_cells;
  double _change = 0.0;

  std::chrono::steady_clock::time_point time_begin;
  std::chrono::steady_clock::time_point time_stamp;
};