// C++ headers
#include <cassert>
#include <iostream>

// Deal.II headers
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

// Project headers
#include <Miscellaneous.h>
#include <Parameter.h>
#include <TopOpt.h>

using namespace dealii;

int main(int argc, char *argv[]) {

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

#ifdef DEBUG
  pcout() << "DEBUG MODE";
#else
  pcout() << "RELEASE MODE";
#endif

  pcout() << " | " << LA_NAME;

  pcout() << " | " << n_mpi_processes() << " PROC";

  std::string parameter_file;

  ParameterManager::getInstance().DeclareParameters();

  if (argc < 2) {
    if (this_mpi_process() == 0)
      ParameterManager::getInstance().OutputDefaultParameters();
    return 0;
  } else if (argc == 2)
    parameter_file = argv[1];
  else
    throw std::runtime_error("Wrong number of parameters");

  ParameterManager::getInstance().ParseParameters(parameter_file);
  auto params = parameters();

  pcout() << " | " << params.general.bvp_type << " " << params.general.problem_type << std::endl;

  TopOpt topopt;
  topopt.Run();

  return 0;

} // End of main
