// C++ headers

// Deal.II headers

// Project headers
#include <Parameter.h>

// DeclareParameters
void ParameterManager::DeclareParameters() {

  // general
  _prm.declare_entry("Analysis name", "topopt", Patterns::Anything(),
                     "Name of the analysis to prefix the output files.");
  _prm.declare_entry("Destination path", "./results/", Patterns::DirectoryName(),
                     "Destination path to write the output files.");
  _prm.declare_entry("Problem type", "opt", Patterns::Anything(),
                     "Options: \"std\": just BVP | \"opt\": optimization.");
  _prm.declare_entry("BVP Type", "elasticityLin", Patterns::Anything(),
                     "Options: \"elasticityLin\": linear elastic");
  _prm.declare_entry("Verbose output", "false", Patterns::Bool(),
                     "Additional console output for debugging.");

  // Domain Geometry
  _prm.enter_subsection("Domain");
  {
    _prm.declare_entry("Geometry name", "hyper_rectangle", Patterns::Anything(),
                       "Options: \"abaqus\": abaqus file | \"hyper_rectangle\": deal.II hyper "
                       "rectangle | \"hyper_L\": deal.II hyper L-beam.");
    _prm.declare_entry("Global refinements", "0", Patterns::Integer(0, 10),
                       "Number of global refinements after mesh generation.");
    _prm.declare_entry(
        "Padding", "false", Patterns::Bool(),
        "Whether to add padding to the domain. If true, the domain will be "
        "extended by 1 cell in each direction to account for filtering boundary effects.");

    _prm.enter_subsection("Abaqus");
    {
      _prm.declare_entry("Input file", "./abq.inp", Patterns::Anything(),
                         "Path to abaqus input file.");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Hyper rectangle");
    {
      _prm.declare_entry("Hyper rectangle point1", "-10,-5", Patterns::List(Patterns::Double()),
                         "Coordinates of point 1: x,y (hyper_rectangle).");
      _prm.declare_entry("Hyper rectangle point2", "10,5", Patterns::List(Patterns::Double()),
                         "Coordinates of point 2: x,y (hyper_rectangle).");
      _prm.declare_entry("Hyper rectangle subdivisions", "20,10",
                         Patterns::List(Patterns::Integer()),
                         "Coarsest level discretization - number of elements per direction: "
                         "x,y,(z) (hyper_rectangle).");
    }
    _prm.leave_subsection(); // leaving the Hyper rectangle subsection.

    _prm.enter_subsection("Hyper L");
    {
      _prm.declare_entry("Hyper L point1", "-10,-10", Patterns::List(Patterns::Double()),
                         "Coordinates of point 1: x,y (hyper_L).");
      _prm.declare_entry("Hyper L point2", "10,10", Patterns::List(Patterns::Double()),
                         "Coordinates of point 2: x,y (hyper_L).");
      _prm.declare_entry("Hyper L subdivisions", "20,20", Patterns::List(Patterns::Integer()),
                         "Subdivision in : x,y (hyper_L).");
      _prm.declare_entry("Hyper L cells to remove", "12,12", Patterns::List(Patterns::Integer()),
                         "Cells to cut off from upper right corner : x,y (hyper_L).");
    }
    _prm.leave_subsection(); // leaving the Hyper L subsection.

    _prm.enter_subsection("Boundary");
    {
      _prm.declare_entry("Boundary IDs", "0,1,1", Patterns::List(Patterns::Integer(0, 50)),
                         "List of boundary IDs to create for each dof, if two boundaries "
                         "satisfy the condition only the first shall survive.");
      _prm.declare_entry("Boundary comp", "0,0,1", Patterns::List(Patterns::Integer(0, 1)),
                         "List of boundary components to be considered while "
                         "deciding the boundary ID.");
      _prm.declare_entry("Boundary coord", "-10,10,0", Patterns::List(Patterns::Double()),
                         "List of coordinates of the boundary, each value "
                         "correspond to the component prescribed above.");
      _prm.declare_entry("Boundary tol", "1e-4,1e-4,1", Patterns::List(Patterns::Double()),
                         "List of tolerance for each boundary coord to satisfy.");
    }
    _prm.leave_subsection(); // leaving boundary subsection.
  }
  _prm.leave_subsection(); // leaving the Domain subsection

  _prm.enter_subsection("BVP");
  {
    _prm.declare_entry("Polynomial degree", "1", Patterns::Integer(1, 3),
                       "Basis function polynomial degree.");

    _prm.enter_subsection("Dirichlet BC");
    {
      _prm.declare_entry("Dirichlet ID", "0,0", Patterns::List(Patterns::Integer(0, 50)),
                         "List of boundary IDs to be assigned Dirichlet BC.");
      _prm.declare_entry("Dirichlet comp", "0,1", Patterns::List(Patterns::Integer(0, 1)),
                         "Components to be assigned Dirichlet BC for respective Dirichlet IDs.");
      _prm.declare_entry("Dirichlet value", "0,0", Patterns::List(Patterns::Double()),
                         "Dirichlet BC values for the respective IDs and components.");
    }
    _prm.leave_subsection(); // leaving Dirichlet BC subsection.

    _prm.enter_subsection("Neumann BC");
    {
      _prm.declare_entry("Neumann ID", "1", Patterns::List(Patterns::Integer(0, 50)),
                         "List of boundary IDs to be assigned Neumann BC.");
      _prm.declare_entry("Neumann comp", "0", Patterns::List(Patterns::Integer(0, 1)),
                         "Components to be assigned Neumann BC for the respective Neumann IDs.");
      _prm.declare_entry("Neumann value", "0", Patterns::List(Patterns::Double()),
                         "Neumann BC values for the respective IDs and components.");
    }
    _prm.leave_subsection(); // leaving Neumann BC subsection.

    _prm.enter_subsection("Newton-Raphson method");
    {
      _prm.declare_entry("Max NR iterations", "100", Patterns::Integer(1, 999),
                         "Max Newton-Raphson iterations.");
      _prm.declare_entry("Tolerance NR convergence", "1e-6", Patterns::Double(),
                         "Tolerance for error in Newton-Raphson method.");
      _prm.declare_entry("Number of load steps", "1", Patterns::Integer(1, 999),
                         "Number of load steps.");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Energy interpolation");
    {
      _prm.declare_entry("Active", "false", Patterns::Bool(),
                         "Use nonlinear-linear energy interpolation in geometrically nonlinear "
                         "topology optimization?");
      _prm.declare_entry("Heaviside beta", "500", Patterns::Double(1, 10000),
                         "Parameter controling the sharpness of heaviside "
                         "projection: def = 500.");
      _prm.declare_entry("Heaviside threshold", "0.01", Patterns::Double(0.0001, 1.0),
                         "Threshold for the heaviside projection: def = 0.01.");
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection(); // leaving BVP subsection.

  _prm.enter_subsection("Material");
  {
    _prm.declare_entry("Material law", "NeoHooke", Patterns::Anything(),
                       "Options: \"StVenant\" | \"NeoHooke\" (default) | \"NeoHookeAD\"");
    _prm.enter_subsection("Elastic isotropic");
    {
      _prm.declare_entry("lambda", "0.5769", Patterns::Double(), "Lamé's first parameter.");
      _prm.declare_entry("mu", "0.3846", Patterns::Double(),
                         "Lamé's second parameter, shear modulus sometimes denoted as G.");
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection(); // leaving material section

  _prm.enter_subsection("Solver");
  {
    _prm.declare_entry("Solver type", "CG", Patterns::Anything(), "Options: \"CG\" | \"Direct\"");
    _prm.declare_entry("Solver tolerance", "1e-9", Patterns::Double(1e-15, 1e15),
                       "Tolerance for iterative solver: actual_tol = "
                       "_system_rhs.l2_norm() * solver_tol.");
  }
  _prm.leave_subsection(); // Leaving solver section.

  _prm.enter_subsection("Optimization");
  {
    _prm.declare_entry("Optimizer", "MMA", Patterns::Anything(), "Options: \"MMA\" (default)");
    _prm.declare_entry("Continuation interval", "10", Patterns::Integer(1, 50),
                       "Number of iterations between each continuation update for all parameters "
                       "that have continuation.");
    _prm.declare_entry("Continuation start iteration", "1", Patterns::Integer(1, 500),
                       "Iteration to start the continuation updates.");
    _prm.declare_entry("Non-design boundary IDs", "", Patterns::List(Patterns::Integer(0, 100)),
                       "Boundary IDs marked as non-design space.");
    _prm.declare_entry("Non-design boundary thickness", "0.0",
                       Patterns::List(Patterns::Double(0, 1e20)),
                       "Thickness of the non-design regions marked by selected boundary IDs.");

    _prm.enter_subsection("GOCM");
    {
      _prm.declare_entry("Step length", "0.05", Patterns::Double(1e-10, 1), "Initial step length.");
      _prm.declare_entry("Step length decay", "0.98", Patterns::Double(0.5, 1),
                         "Multiplier for the step length decay for every iteration.");
      _prm.declare_entry("Minimum step length", "0.001", Patterns::Double(0, 1),
                         "Minimum possible step length.");
      _prm.declare_entry("Update formula", "exponential", Patterns::Anything(),
                         "Options: \"exponential\" (default) | \"linear\" | \"reciprocal\"");
      _prm.declare_entry("Move limit lagrange", "0.5", Patterns::Double(0.05, 1),
                         "Multiplicative move limit for the lagrange multipliers in the form: "
                         "Lam^(i+1) = min(Lam^i * (1 - move), max(Lam^i * (1 + move), Lam^(i+1)))");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("MMA");
    {
      _prm.declare_entry("Move limit", "0.1", Patterns::Double(1e-10, 1e20), "MMA move limit");
      _prm.declare_entry("Asymptote move initial", "0.1", Patterns::Double(1e-10, 1e20),
                         "Asymptote move initial.");
      _prm.declare_entry("Asymptote move decrease", "0.7", Patterns::Double(1e-10, 1e20),
                         "Asymptote move decrease.");
      _prm.declare_entry("Asymptote move increase", "1.2", Patterns::Double(1e-10, 1e20),
                         "Asymptote move increase.");
      _prm.declare_entry("Robust asymptotes type", "0", Patterns::Integer(0, 1),
                         "Robust asymptotes type: [0, 1]");
      _prm.declare_entry(
          "Constraint modification", "false", Patterns::Bool(),
          "Alternative way to compute pij and qij in MMA. Options: \"false\" (default) | \"true\"");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Design field");
    {
      _prm.declare_entry(
          "Mode", "simp", Patterns::Anything(),
          "Options: \"simp\" (default) | \"variable-thickness\" | \"multi-thickness\"");
      _prm.declare_entry("Initial density", "0.5", Patterns::Double(0.001, 1.0),
                         "Initial value of densities.");
      _prm.declare_entry("Penalty", "3.0", Patterns::Double(1.0, 100.0),
                         "Penalization exponent for intermediate densities.");
      _prm.declare_entry("Stress relaxation epsilon initial", "0.1", Patterns::Double(1e-20, 1.0),
                         "Epsilon-relaxation parameter for the stresses.");
      _prm.declare_entry("Stress relaxation epsilon target", "0.01", Patterns::Double(1e-20, 1.0),
                         "Target value for the epsilon-relaxation parameter for the stresses.");
      _prm.declare_entry("Stress relaxation target iteration", "150", Patterns::Integer(1.0, 500.0),
                         "Target iteration to reach the final epsilon-relaxation error.");
      _prm.enter_subsection("Variable-thickness");
      {
        _prm.declare_entry("Penalty continuation", "false", Patterns::Bool(),
                           "Perform penalty continuation from 1 to penalty value?");
        _prm.declare_entry(
            "Penalty target iteration", "50", Patterns::Integer(1.0, 500.0),
            "Target iteration for the penalization exponent for intermediate densities");
        _prm.declare_entry("Low thickness threshold", "0.1", Patterns::Double(0.0, 1.0),
                           "Thin features below this density value will be eliminated.");
      }
      _prm.leave_subsection();
      _prm.enter_subsection("Multi-thickness");
      {
        _prm.declare_entry(
            "Penalty target iteration", "100", Patterns::Integer(1, 500),
            "Target iteration to reach the final penalization exponent for intermediate densities");
        _prm.declare_entry("Number of thicknesses", "2", Patterns::Integer(1, 10),
                           "Number of target thicknesses.");
      }
      _prm.leave_subsection();
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Filters");
    {
      _prm.declare_entry("Blur filter", "true", Patterns::Bool(), "Apply density filter?");
      _prm.declare_entry("Blur filter type", "radius", Patterns::Anything(),
                         "Options: \"PDE\" (default) | \"radius\"");
      _prm.declare_entry("Blur filter radius", "0.2", Patterns::Double(1e-20, 1e6),
                         "Value of the density filter radius.");
      _prm.declare_entry(
          "PDE filter boundary penalization", "false", Patterns::Bool(),
          "Apply boundary penalization in PDE filter to eliminate boundary \"sticking\"?");
      _prm.declare_entry("Heaviside projection", "true", Patterns::Bool(),
                         "Apply Heaviside projection?");
      _prm.declare_entry("Heaviside projection threshold", "0.5", Patterns::Double(0.001, 1),
                         "Heaviside projection threshold.");
      _prm.declare_entry("Heaviside projection beta", "50", Patterns::Double(1, 1000),
                         "Heaviside projection beta.");
      _prm.declare_entry("Heaviside projection target iteration", "100", Patterns::Integer(1, 500),
                         "Target iteration to reach the final Heaviside projection beta.");
      _prm.enter_subsection("Variable-thickness");
      {
        _prm.declare_entry("DGI projection", "true", Patterns::Bool(), "Apply DGI projection?");
        _prm.declare_entry("DGI projection beta", "10", Patterns::Double(0.1, 1000),
                           "Projection sharpness for the DGI deblurring after filtering");
        _prm.declare_entry("DGI projection target iteration", "100", Patterns::Integer(1, 500),
                           "Target iteration to reach the final DGI projection beta.");
      }
      _prm.leave_subsection();
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Convergence");
    {
      _prm.declare_entry("Design change tolerance", "0.001", Patterns::Double(1e-20, 1),
                         "Design change tolerance.");
      _prm.declare_entry("Max iterations", "100", Patterns::Integer(1, 5000),
                         "Max number of max iterations [1, 5000]");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Adaptive meshing");
    {
      _prm.declare_entry("Active", "false", Patterns::Bool(), "Use adaptive meshing?");
      _prm.declare_entry("Start iteration", "5", Patterns::Integer(0, 500),
                         "At which iteration should the first mesh adaptation take place?");
      _prm.declare_entry("Adaptation interval", "5", Patterns::Integer(1, 500),
                         "What is iteration interval for the mesh adaptation?");
      _prm.declare_entry("Criterion", "DENSITY-JUMP", Patterns::Anything(),
                         "DENSITY-JUMP | DENSITY | CNF | VONMISES");
      _prm.declare_entry("Max h refinements", "3", Patterns::Integer(0, 10),
                         "Max h refinement level allowed");
      _prm.declare_entry("Min h refinements", "0", Patterns::Integer(0, 10),
                         "Min h refinement level allowed");
      _prm.declare_entry("Relative threshold refinement", "0.25", Patterns::Double(1.0e-20, 1.0),
                         "Relative threshold is a multiplier for maximum criterion value.");
      _prm.declare_entry("Relative threshold coarsening", "1e-2", Patterns::Double(1.0e-20, 1e-1),
                         "Relative threshold is a multiplier for maximum criterion value.");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Responses");
    {
      _prm.declare_entry(
          "Objective", "compliance", Patterns::Anything(),
          "Options (minimization): \"compliance\" (default) | \"volume\" | \"PNStress\".");
      _prm.declare_entry("Constraints", "volume", Patterns::List(Patterns::Anything()),
                         "Options (list): \"volume\" (default) | \"compliance\" | \"PNStress\".");
      _prm.declare_entry("Constraints types", "upper", Patterns::List(Patterns::Anything()),
                         "Options (list): \"upper\", \"lower\", \"equality\".");
      _prm.declare_entry("Constraints values", "0.5", Patterns::List(Patterns::Double(-1e20, 1e20)),
                         "Bounds for each constraint (list).");
      _prm.enter_subsection("Von-Mises options");
      {
        _prm.declare_entry("Allowable Von-Mises", "1", Patterns::Double(1e-20, 1e20),
                           "Allowable Von-Mises stress.");
        _prm.declare_entry("P-Norm exponent", "10", Patterns::Double(1, 1e20),
                           "Initial value of P in P-norm.");
      }
      _prm.leave_subsection();
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection(); // Leaving optimization section

} // End of DeclareParameters function

void ParameterManager::ParseParameters(std::string filename) {

  _prm.parse_input(filename);
  std::string arg;

  _parameters.file.name = _prm.get("Analysis name");
  _parameters.file.destination_path = _prm.get("Destination path");
  _parameters.general.problem_type = _prm.get("Problem type");
  _parameters.general.bvp_type = _prm.get("BVP Type");

  _parameters.debug.verbose = _prm.get_bool("Verbose output");

  _prm.enter_subsection("Domain");
  {
    _parameters.geometry.name = _prm.get("Geometry name");
    _parameters.geometry.global_refinements = _prm.get_integer("Global refinements");
    _parameters.geometry.padding = _prm.get_bool("Padding");

    _prm.enter_subsection("Abaqus");
    {
      _parameters.geometry.abaqus.input = _prm.get("Input file");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Hyper rectangle");
    {
      arg = _prm.get("Hyper rectangle point1");
      string_to_vector_of_double(arg, _parameters.geometry.hyper_rectangle.point1);

      arg = _prm.get("Hyper rectangle point2");
      string_to_vector_of_double(arg, _parameters.geometry.hyper_rectangle.point2);

      arg = _prm.get("Hyper rectangle subdivisions");
      string_to_vector_of_int(arg, _parameters.geometry.hyper_rectangle.subdivision);
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Hyper L");
    {
      arg = _prm.get("Hyper L point1");
      string_to_vector_of_double(arg, _parameters.geometry.hyper_l.point1);

      arg = _prm.get("Hyper L point2");
      string_to_vector_of_double(arg, _parameters.geometry.hyper_l.point2);

      arg = _prm.get("Hyper L subdivisions");
      string_to_vector_of_int(arg, _parameters.geometry.hyper_l.subdivision);

      arg = _prm.get("Hyper L cells to remove");
      string_to_vector_of_int(arg, _parameters.geometry.hyper_l.cells_to_remove);
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Boundary");
    {
      arg = _prm.get("Boundary IDs");
      string_to_vector_of_int(arg, _parameters.geometry.boundary_ids);

      arg = _prm.get("Boundary comp");
      string_to_vector_of_int(arg, _parameters.geometry.boundary_comp);

      arg = _prm.get("Boundary coord");
      string_to_vector_of_double(arg, _parameters.geometry.boundary_coord);

      arg = _prm.get("Boundary tol");
      string_to_vector_of_double(arg, _parameters.geometry.boundary_tol);
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection();

  _prm.enter_subsection("BVP");
  {
    _parameters.bvp.poly_degree = _prm.get_integer("Polynomial degree");

    // dirichlet_bc
    _prm.enter_subsection("Dirichlet BC");
    {
      arg = _prm.get("Dirichlet ID");
      string_to_vector_of_int(arg, _parameters.bvp.dirichlet.id);

      arg = _prm.get("Dirichlet comp");
      string_to_vector_of_int(arg, _parameters.bvp.dirichlet.component);

      arg = _prm.get("Dirichlet value");
      string_to_vector_of_double(arg, _parameters.bvp.dirichlet.value);

      if (_parameters.bvp.dirichlet.id.size() != _parameters.bvp.dirichlet.component.size() ||
          _parameters.bvp.dirichlet.id.size() != _parameters.bvp.dirichlet.value.size())
        throw std::invalid_argument("Unequal number of parameters for Dirichlet BC");
    }
    _prm.leave_subsection(); // leaving Dirichlet BC subsection.

    // neumann_bc
    _prm.enter_subsection("Neumann BC");
    {
      arg = _prm.get("Neumann ID");
      string_to_vector_of_int(arg, _parameters.bvp.neumann.id);

      arg = _prm.get("Neumann comp");
      string_to_vector_of_int(arg, _parameters.bvp.neumann.component);

      arg = _prm.get("Neumann value");
      string_to_vector_of_double(arg, _parameters.bvp.neumann.value);

      if (_parameters.bvp.neumann.id.size() != _parameters.bvp.neumann.component.size() ||
          _parameters.bvp.neumann.id.size() != _parameters.bvp.neumann.value.size())
        throw std::invalid_argument("Unequal number of parameters for Neumann BC");
    }
    _prm.leave_subsection(); // leaving Neumann BC subsection.

    // Newton-Raphson method
    _prm.enter_subsection("Newton-Raphson method");
    {
      _parameters.bvp.nonlinear.max_NR_iter = _prm.get_integer("Max NR iterations");
      _parameters.bvp.nonlinear.tol_residual_NR = _prm.get_double("Tolerance NR convergence");
      _parameters.bvp.nonlinear.number_load_steps = _prm.get_integer("Number of load steps");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Energy interpolation");
    {
      _parameters.bvp.nonlinear.energy_interpolation.active = _prm.get_bool("Active");
      _parameters.bvp.nonlinear.energy_interpolation.heaviside_beta =
          _prm.get_double("Heaviside beta");
      _parameters.bvp.nonlinear.energy_interpolation.heaviside_threshold =
          _prm.get_double("Heaviside threshold");
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection(); // leaving BVP subsection.

  // Material
  _prm.enter_subsection("Material");
  {
    _parameters.material.law = _prm.get("Material law");
    _prm.enter_subsection("Elastic isotropic");
    {
      _parameters.material.elastic_isotropic.lambda = _prm.get_double("lambda");
      _parameters.material.elastic_isotropic.mu = _prm.get_double("mu");
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection(); // leaving material section

  // Solver
  _prm.enter_subsection("Solver");
  {
    _parameters.bvp.solver.name = _prm.get("Solver type");
    _parameters.bvp.solver.tol = _prm.get_double("Solver tolerance");
  }
  _prm.leave_subsection(); // leaving solver section

  // Optimization
  _prm.enter_subsection("Optimization");
  {
    _parameters.optimization.optimizer = _prm.get("Optimizer");
    _parameters.optimization.continuation_interval = _prm.get_integer("Continuation interval");
    _parameters.optimization.continuation_start_iteration =
        _prm.get_integer("Continuation start iteration");
    arg = _prm.get("Non-design boundary IDs");
    string_to_vector_of_int(arg, _parameters.optimization.non_design_id);
    arg = _prm.get("Non-design boundary thickness");
    string_to_vector_of_double(arg, _parameters.optimization.non_design_thickness);
    _prm.enter_subsection("GOCM");
    {
      _parameters.optimization.gocm.step_length = _prm.get_double("Step length");
      _parameters.optimization.gocm.step_length_decay = _prm.get_double("Step length decay");
      _parameters.optimization.gocm.min_step_length = _prm.get_double("Minimum step length");
      _parameters.optimization.gocm.update_formula = _prm.get("Update formula");
      _parameters.optimization.gocm.move_limit_lagrange = _prm.get_double("Move limit lagrange");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("MMA");
    {
      _parameters.optimization.mma.move_limit = _prm.get_double("Move limit");
      _parameters.optimization.mma.asym_move_init = _prm.get_double("Asymptote move initial");
      _parameters.optimization.mma.asym_move_dec = _prm.get_double("Asymptote move decrease");
      _parameters.optimization.mma.asym_move_incr = _prm.get_double("Asymptote move increase");
      _parameters.optimization.mma.robust_asymptotes_type =
          _prm.get_integer("Robust asymptotes type");
      _parameters.optimization.mma.constraint_modification =
          _prm.get_bool("Constraint modification");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Design field");
    {
      _parameters.optimization.design_field.mode = _prm.get("Mode");
      _parameters.optimization.design_field.initial_density = _prm.get_double("Initial density");
      _parameters.optimization.design_field.penalty = _prm.get_double("Penalty");
      _parameters.optimization.design_field.stress_relaxation_epsilon_initial =
          _prm.get_double("Stress relaxation epsilon initial");
      _parameters.optimization.design_field.stress_relaxation_epsilon_target =
          _prm.get_double("Stress relaxation epsilon target");
      _parameters.optimization.design_field.stress_relaxation_target_iteration =
          _prm.get_integer("Stress relaxation target iteration");
      _prm.enter_subsection("Variable-thickness");
      {
        _parameters.optimization.design_field.variable_thickness.penalty_continuation =
            _prm.get_bool("Penalty continuation");
        _parameters.optimization.design_field.variable_thickness.penalty_target_iteration =
            _prm.get_integer("Penalty target iteration");
        _parameters.optimization.design_field.variable_thickness.low_thickness_threshold =
            _prm.get_double("Low thickness threshold");
      }
      _prm.leave_subsection();
      _prm.enter_subsection("Multi-thickness");
      {
        _parameters.optimization.design_field.multi_thickness.penalty_target_iteration =
            _prm.get_integer("Penalty target iteration");
        _parameters.optimization.design_field.multi_thickness.n_thicknesses =
            _prm.get_double("Number of thicknesses");
      }
      _prm.leave_subsection();
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Filters");
    {
      _parameters.optimization.filters.blur_filter = _prm.get_bool("Blur filter");
      _parameters.optimization.filters.blur_filter_type = _prm.get("Blur filter type");
      _parameters.optimization.filters.blur_filter_radius = _prm.get_double("Blur filter radius");
      _parameters.optimization.filters.pde_boundary_penalization =
          _prm.get_bool("PDE filter boundary penalization");
      _parameters.optimization.filters.heaviside_projection = _prm.get_bool("Heaviside projection");
      _parameters.optimization.filters.heaviside_projection_threshold =
          _prm.get_double("Heaviside projection threshold");
      _parameters.optimization.filters.heaviside_projection_beta =
          _prm.get_double("Heaviside projection beta");
      _parameters.optimization.filters.heaviside_projection_target_iteration =
          _prm.get_integer("Heaviside projection target iteration");
      _prm.enter_subsection("Variable-thickness");
      {
        _parameters.optimization.filters.variable_thickness.dgi_projection =
            _prm.get_bool("DGI projection");
        _parameters.optimization.filters.variable_thickness.dgi_projection_beta =
            _prm.get_double("DGI projection beta");
        _parameters.optimization.filters.variable_thickness.dgi_projection_target_iteration =
            _prm.get_integer("DGI projection target iteration");
      }
      _prm.leave_subsection();
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Convergence");
    {
      _parameters.optimization.convergence.design_change_tol =
          _prm.get_double("Design change tolerance");
      _parameters.optimization.convergence.max_iterations = _prm.get_integer("Max iterations");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Adaptive meshing");
    {
      _parameters.optimization.mesh_adaptivity.active = _prm.get_bool("Active");
      _parameters.optimization.mesh_adaptivity.start_iteration =
          _prm.get_integer("Start iteration");
      _parameters.optimization.mesh_adaptivity.interval_iteration =
          _prm.get_integer("Adaptation interval");
      _parameters.optimization.mesh_adaptivity.type = _prm.get("Criterion");
      _parameters.optimization.mesh_adaptivity.max_h_refinements =
          _prm.get_integer("Max h refinements");
      _parameters.optimization.mesh_adaptivity.min_h_refinements =
          _prm.get_integer("Min h refinements");
      _parameters.optimization.mesh_adaptivity.relative_threshold_refine =
          _prm.get_double("Relative threshold refinement");
      _parameters.optimization.mesh_adaptivity.relative_threshold_coarsen =
          _prm.get_double("Relative threshold coarsening");
    }
    _prm.leave_subsection();

    _prm.enter_subsection("Responses");
    {
      unsigned int response_id = 0;
      auto enum_val = stringToEnum<ResponseNames>(_prm.get("Objective"), ResponseMap);
      if (enum_val.has_value()) {
        _parameters.optimization.objective.name = enum_val.value();
        _parameters.optimization.objective.id = response_id;
        ++response_id;
      } else
        throw std::invalid_argument("Wrong objective name given in the parameter file.");

      std::vector<std::string> vec_str;

      arg = _prm.get("Constraints");
      string_to_vector_of_strings(arg, vec_str);
      _parameters.optimization.constraints.resize(vec_str.size());
      for (unsigned int i = 0; i < vec_str.size(); ++i) {
        auto enum_val = stringToEnum<ResponseNames>(vec_str[i], ResponseMap);
        if (enum_val.has_value()) {
          _parameters.optimization.constraints[i].name = enum_val.value();
          _parameters.optimization.constraints[i].id = response_id;
          ++response_id;
        } else
          throw std::invalid_argument("Wrong constraint name given in the parameter file.");
      }

      arg = _prm.get("Constraints types");
      string_to_vector_of_strings(arg, vec_str);
      if (vec_str.size() == _parameters.optimization.constraints.size()) {
        for (unsigned int i = 0; i < vec_str.size(); ++i) {
          auto enum_val = stringToEnum<ConstraintTypeNames>(vec_str[i], ConstraintTypeMap);
          if (enum_val.has_value())
            _parameters.optimization.constraints[i].type = enum_val.value();
          else
            throw std::invalid_argument("Wrong response type given in the parameter file.");
        }
      } else
        throw std::invalid_argument("Wrong number of constraint types in the parameter file. ");

      std::vector<double> vec_double;
      arg = _prm.get("Constraints values");
      string_to_vector_of_double(arg, vec_double);
      if (vec_double.size() == _parameters.optimization.constraints.size())
        for (unsigned int i = 0; i < vec_double.size(); ++i)
          _parameters.optimization.constraints[i].constr_value = vec_double[i];
      else
        throw std::invalid_argument("Wrong number of constraint values in the parameter file. ");
    }
    _prm.leave_subsection();
  }
  _prm.leave_subsection();

} // End of ParseParameters fucntion

void ParameterManager::OutputDefaultParameters() {
  const std::string out_file = "default.prm";
  std::ofstream parameter_out(out_file);
  _prm.print_parameters(parameter_out, ParameterHandler::Text);
}