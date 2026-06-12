#pragma once

// C++ headers
#include <assert.h>
#include <fstream>
#include <iostream>

// Deal.II headers
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>

// Project headers
#include <MacrosAndTypedefs.h>
#include <Miscellaneous.h>

using namespace dealii;

struct FileParameters {
  std::string name = "topology-optimization";  /// Any name to prefix the output files
  std::string destination_path = "./results/"; /// Path to write the output files
};

struct GeneralParameters {
  std::string problem_type = "opt";             /// bvp | opt; i.e standard or optimization
  std::string bvp_type = "NonlinearElasticity"; /// LinearElasticity | NonlinearElasticity
  bool use_automatic_differentiation = false;   /// use automatic differentiation
};

struct DebugParameters {
  bool verbose = false; // Additional output for debugging. false == no output;
                        // true == yes.
};

struct AbaqusParameters {
  std::string input;
};

struct HyperRectangleParameters {
  std::vector<double> point1, // bottom left corner
      point2;                 // top right corner
  std::vector<int> subdivision = {1, 1, 1};
};

struct HyperLParameters {
  std::vector<double> point1, // bottom left corner
      point2;                 // top right corner
  std::vector<int> subdivision = {1, 1, 1};
  std::vector<int> cells_to_remove = {1, 1, 1};
};

struct GeometryParameters {
  std::string name = "hyper_rectangle";       // abaqus | hyper_rectangle | hyper_cube_with_hole
  std::vector<int> boundary_ids = {1, 2, 2};  // ID to assign to a bounary
  std::vector<int> boundary_comp = {0, 0, 1}; // components to locate the boundary
  std::vector<double> boundary_coord = {-1, 1, 0}; // coordinates to locate the boundary
  std::vector<double> boundary_tol = {
      1e-3, 1e-3, 0.1}; // tolerance within which a point should be in relation to the
                        // boundary_coord to be assigned a particular boundary ID.

  AbaqusParameters abaqus;
  HyperRectangleParameters hyper_rectangle;
  HyperLParameters hyper_l;
  unsigned int global_refinements = 0;
  bool padding = false;
};

struct BoundaryConditionsParameters {
  std::vector<int> id;
  std::vector<int> component;
  std::vector<double> value;
};

struct EnergyInterpolationParameters {
  bool active = false;
  double heaviside_beta = 500.0;
  double heaviside_threshold = 0.01;
};

struct SolverParameters {
  std::string name = "CG";
  int max_iter = 100;
  double tol = 1e-4;
};

struct NonlinearBVPParameters {
  int max_NR_iter = 1;
  double tol_residual_NR = 1e-8;
  int number_load_steps = 1;
  EnergyInterpolationParameters energy_interpolation;
};

struct BVPParameters {
  unsigned int poly_degree = 1;
  BoundaryConditionsParameters dirichlet;
  BoundaryConditionsParameters neumann;
  SolverParameters solver;
  NonlinearBVPParameters nonlinear;
};

struct ElasticIsotropicParameters {
  double lambda = 115000.0;
  double mu = 77000;
};

struct MaterialParameters {
  std::string law = "NeoHooke";
  ElasticIsotropicParameters elastic_isotropic;
};

struct VariableThicknessParameters {
  bool penalty_continuation = false;
  unsigned int penalty_target_iteration = 50;
  double low_thickness_threshold = 0.05;
};

struct MultiThicknessParameters {
  unsigned int penalty_target_iteration = 50;
  unsigned int n_thicknesses = 2;
};

struct DesignFieldParameters {
  std::string mode = "simp";
  double initial_density = 0.5;
  double penalty = 3.0;
  double stress_relaxation_epsilon_initial = 0.1;
  double stress_relaxation_epsilon_target = 0.01;
  unsigned int stress_relaxation_target_iteration = 100;
  VariableThicknessParameters variable_thickness;
  MultiThicknessParameters multi_thickness;
};

struct FilterVariableThicknessParameters {
  bool dgi_projection = true;
  double dgi_projection_beta = 25;
  unsigned int dgi_projection_target_iteration = 100;
};

struct FilterParameters {
  bool geometry_constraints = false;
  std::vector<std::string> geometry_constraints_types;
  bool heaviside_projection = true;
  double heaviside_projection_threshold = 0.5;
  double heaviside_projection_beta = 25;
  unsigned int heaviside_projection_target_iteration = 100;
  bool blur_filter = true;
  std::string blur_filter_type = "PDE";
  double blur_filter_radius = 0.5;
  bool pde_boundary_penalization = false;
  FilterVariableThicknessParameters variable_thickness;
};

struct MMAParameters {
  double move_limit = 0.2;
  double asym_move_init = 0.5;
  double asym_move_dec = 0.7;
  double asym_move_incr = 1.2;
  unsigned int robust_asymptotes_type = 0;
  bool constraint_modification = false;
};

struct GOCMParameters {
  double step_length = 0.05;
  double step_length_decay = 0.98;
  double min_step_length = 0.001;
  std::string update_formula;
  double move_limit_lagrange = 0.5;
};

struct ResponseParameters {
  ResponseNames name;
  unsigned int id;
};

struct ConstraintParameters : public ResponseParameters {
  ConstraintTypeNames type = ConstraintTypeNames::upper;
  double constr_value = 1.0;
};

struct MeshAdaptivityParameters {
  bool active = false;
  unsigned int start_iteration;
  unsigned int interval_iteration;
  std::string type;
  unsigned int max_h_refinements;
  unsigned int min_h_refinements;
  double relative_threshold_refine;
  double relative_threshold_coarsen;
};

struct ConvergenceParameters {
  double design_change_tol = 0.0;
  unsigned int max_iterations;
};

struct OptimizationParameters {
  std::string optimizer;
  unsigned int continuation_interval = 1;
  unsigned int continuation_start_iteration = 1;
  std::vector<int> non_design_id;
  std::vector<double> non_design_thickness;
  DesignFieldParameters design_field;
  FilterParameters filters;
  MMAParameters mma;
  GOCMParameters gocm;
  ResponseParameters objective;
  std::vector<ConstraintParameters> constraints;
  ConvergenceParameters convergence;
  MeshAdaptivityParameters mesh_adaptivity;
};

struct Parameters {
  FileParameters file;
  GeneralParameters general;
  DebugParameters debug;
  GeometryParameters geometry;
  BVPParameters bvp;
  MaterialParameters material;
  OptimizationParameters optimization;
};

class ParameterManager {
public:
  static ParameterManager &getInstance() {
    static ParameterManager instance;
    return instance;
  }

  const Parameters &getParameters() { return _parameters; };
  void DeclareParameters();
  void ParseParameters(std::string filename = "stop.prm");
  void OutputDefaultParameters();

private:
  ParameterManager() {} // Private constructor for singleton
  ParameterManager(ParameterManager const &) = delete;
  void operator=(ParameterManager const &) = delete;

  Parameters _parameters;
  ParameterHandler _prm;
};

inline const Parameters &parameters() { return ParameterManager::getInstance().getParameters(); }