#ifndef MACROS_AND_TYPEDEFS_H
#define MACROS_AND_TYPEDEFS_H

// C++ headers
#include <fstream>
#include <iostream>

// Deal.II headers
#include <deal.II/algorithms/general_data_storage.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/differentiation/ad.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

// Project headers
#include "GlobalContext.h"

// Your macros, type aliases, and other definitions here.
namespace LA {
#ifdef USE_PETSC
using namespace dealii::LinearAlgebraPETSc;
#define LA_NAME "PETSc"
#else
using namespace dealii::LinearAlgebraTrilinos;
#define LA_NAME "Trilinos"
#endif
} // namespace LA

/// Macro to have an MPI_Barrier
#define PROJ_MPI_BARRIER                                                                           \
  {                                                                                                \
    int ierr = MPI_Barrier(mpi_communicator());                                                    \
    AssertThrowMPI(ierr);                                                                          \
  }

#define FILE_LINE                                                                                  \
  {                                                                                                \
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;                                         \
  }

#define FILE_LINE_P                                                                                \
  {                                                                                                \
    std::cout << __FILE__ << ":" << __LINE__ << " |P: " << this_mpi_process() << std::endl;        \
  }

using namespace dealii;

using DomainParallelTriaType = parallel::TriangulationBase<2>; // base background mesh

using DomainParallelTriaTypeDistributed = parallel::distributed::Triangulation<2>;

struct DataComponentInterpretationTypes {
  static inline std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
      scalar_interpretation = {dealii::DataComponentInterpretation::component_is_scalar};

  static inline std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
      vec_interpretation = {
          std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>(
              2, dealii::DataComponentInterpretation::component_is_part_of_vector)};
};

class StandardTensors {
public:
  static const SymmetricTensor<2, 2> I() { return unit_symmetric_tensor<2>(); }
  static const SymmetricTensor<4, 2> IxI() { return outer_product(I(), I()); }
  static const SymmetricTensor<4, 2> II() { return identity_tensor<2>(); }
  static const SymmetricTensor<4, 2> II_dev() { return deviator_tensor<2>(); }
};

constexpr Differentiation::AD::NumberTypes ADTypeCodeScalar =
    Differentiation::AD::NumberTypes::sacado_dfad_dfad;

using ADScalarHelper =
    Differentiation::AD::ScalarFunction<2, Differentiation::AD::NumberTypes::sacado_dfad_dfad,
                                        double>;
using ADScalarNumberType = typename ADScalarHelper::ad_type;

using ADEnergyHelper =
    Differentiation::AD::EnergyFunctional<Differentiation::AD::NumberTypes::sacado_dfad_dfad,
                                          double>;
using ADEnergyNumberType = typename ADEnergyHelper::ad_type;

namespace OptModeNames {
const std::string VariableThickness = "variable-thickness";
} // namespace OptModeNames

namespace FilterNames {
const std::string PDEFilter = "PDEFilter";
const std::string HeavisideProjection = "HeavisideProjection";
const std::string DGIProjection = "DGIProjection";
const std::string StressFilter = "StressFilter";
const std::string UnknownFilter = "UnknownFilter";
} // namespace FilterNames

namespace DesignFieldNames {
const std::string PseudoDensity = "PseudoDensity";
const std::string FilteredDensity = "FilteredDensity";
const std::string DGIProjectedDensity = "DGIProjectedDensity";
const std::string HeavisideProjectedDensity = "HeavisideProjectedDensity";
const std::string PreviousPseudoDensity = "PreviousPseudoDensity";
} // namespace DesignFieldNames

namespace MMATransferDataNames {
const std::string OldPseudoDensity = "OldPseudoDensity";
const std::string Old2PseudoDensity = "Old2PseudoDensity";
const std::string MinDensity = "MinDensity";
const std::string MaxDensity = "MaxDensity";
const std::string LowerAsymptotes = "LowerAsymptotes";
const std::string UpperAsymptotes = "UpperAsymptotes";
} // namespace MMATransferDataNames

enum class ResponseNames { compliance, volume, PNStress };

static const std::map<std::string, ResponseNames> ResponseMap = {
    {"compliance", ResponseNames::compliance},
    {"volume", ResponseNames::volume},
    {"PNStress", ResponseNames::PNStress}};

enum class ConstraintTypeNames { equality, lower, upper };

static const std::map<std::string, ConstraintTypeNames> ConstraintTypeMap = {
    {"equality", ConstraintTypeNames::equality},
    {"lower", ConstraintTypeNames::lower},
    {"upper", ConstraintTypeNames::upper}};

// Define other enum for postprocessing of other data, e.g. e_cauchy_11.
enum { e_von_mises_s = 100 };

enum { e_padding_mat_id = 999 };

#endif // MACROS_AND_TYPEDEFS_H