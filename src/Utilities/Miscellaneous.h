#pragma once

#include <fstream>
#include <iostream>
#include <optional>

#include <boost/algorithm/string.hpp> // to trim string

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/algorithms/general_data_storage.h>
#include <deal.II/differentiation/ad.h>

// Project file
#include <GlobalContext.h>

using namespace dealii;

/// string to vector of int
void string_to_vector_of_int(const std::string &s, std::vector<int> &int_vec);

/// string to vector of double
void string_to_vector_of_double(const std::string &s, std::vector<double> &double_vec);

/// string to vector of strings
void string_to_vector_of_strings(const std::string &s, std::vector<std::string> &string_vec);

template <typename K, typename V>
inline K findKeyByValue(const std::map<K, V> &map, const V &value) {
  for (const auto &pair : map) {
    if (pair.second == value) {
      return pair.first;
    }
  }
  // Value not found
  throw std::runtime_error("Value not found in map"); // Or return a optional, or a bool, or a
                                                      // default value.
}

/// Type-safe signum function, returns -1,0,1
template <typename T> inline int sgn(T val) { return (T(0) < val) - (val < T(0)); }

std::pair<double, double> cartesianToPolar(const Tensor<1, 2> &cartesian);

class Errors {
  /**
   * @brief This class is to compute the relative error, which is used in Newton
   * Raphson scheme. The class is Initialized with error in the first iteration,
   * subsequently in later iterations given an error value, it can return
   * normalized error wrt the error in first iteration.
   */
private:
  double error_first_iter = 0.0;
  bool initialized = false;

public:
  /**
   * @brief Initialize error.
   *
   * @param error
   */
  inline void Initialize(double error) {
    if (error == 0.0)
      throw std::runtime_error("First iteration error cannot be 0.0 ");
    else {
      if (!initialized) {
        error_first_iter = error;
        initialized = true;
      } else
        std::cerr << "Already the error is initialized." << std::endl;
    }
  }

  /**
   * @brief Function to get the Normalized Error.
   *
   * @param error
   * @return double
   */
  inline double GetNormalizedError(double error) {
    if (initialized)
      return error / error_first_iter;
    else {
      std::cerr << "First iteration error not initialized, so cannot Normalize." << std::endl;
      return 1e9;
    }
  }

  /**
   * @brief Reset the error.
   *
   */
  inline void Reset() {
    error_first_iter = 0.0;
    initialized = false;
  }
};

Vector<double> GetMPIReducedVectorFromMap(std::map<CellId, double> &map);
double L1NormMap(std::map<CellId, double> &map);
double L1NormSizeDividedMap(std::map<CellId, double> &map);
double L2NormMap(std::map<CellId, double> &map);
double L2NormSizeDividedMap(std::map<CellId, double> &map);
double LInfinityNormMap(std::map<CellId, double> &map);
void L1NormalizeMap(std::map<CellId, double> &map);
void L1NormalizeSizeIndependentMap(std::map<CellId, double> &map);
void L2NormalizeMap(std::map<CellId, double> &map);
void L2NormalizeSizeIndependentMap(std::map<CellId, double> &map);
void LInfinityNormalizeMap(std::map<CellId, double> &map);
double MapDotProduct(const std::map<CellId, double> &map_data1,
                     const std::map<CellId, double> &map_data2);
double SmoothedHeavisideProjection(const double val, const double beta = 500,
                                   const double rho_0 = 0.01);
double SmoothedHeavisideProjectionDer(const double val, const double beta = 500,
                                      const double rho_0 = 0.01);

template <typename EnumType>
std::optional<EnumType> stringToEnum(const std::string &str,
                                     const std::map<std::string, EnumType> &stringMap) {
  if (auto it = stringMap.find(str); it != stringMap.end()) {
    return it->second;
  } else {
    return std::nullopt;
  }
}

double SetLBeamNotchDensity(const Point<2> &pt, const double notch_radius, const bool antialiasing,
                            const double cell_length);