#include <Miscellaneous.h>
#include <Parameter.h>

// string to vector of int
void string_to_vector_of_int(const std::string &s, std::vector<int> &int_vec) {
  int_vec.clear();
  std::string item;
  std::stringstream ss(s);

  while (std::getline(ss, item, ',')) {
    int_vec.push_back(std::stoi(item));
  }
}

// string to vector of double
void string_to_vector_of_double(const std::string &s, std::vector<double> &double_vec) {
  double_vec.clear();
  std::string item;
  std::stringstream ss(s);

  while (std::getline(ss, item, ',')) {
    double_vec.push_back(std::stof(item));
  }
}

// string to vector of double
void string_to_vector_of_strings(const std::string &s, std::vector<std::string> &string_vec) {
  string_vec.clear();
  std::string item;
  std::stringstream ss(s);

  while (std::getline(ss, item, ',')) {
    boost::algorithm::trim(item);
    string_vec.push_back(item);
  }
}

std::pair<double, double> cartesianToPolar(const Tensor<1, 2> &cartesian) {
  double x = cartesian[0];
  double y = cartesian[1];

  double r = std::sqrt(x * x + y * y); // Radius
  double theta = std::atan2(y, x);     // Angle (in radians)

  return {r, theta};
}

Vector<double> GetMPIReducedVectorFromMap(std::map<CellId, double> &map) {
  std::vector<std::map<CellId, double>> vec_map =
      Utilities::MPI::all_gather(mpi_communicator(), map);
  std::map<CellId, double> map_reduced;
  for (const auto &map_proc : vec_map)
    for (const auto &[cellid, val] : map_proc) {
      if (map_reduced.find(cellid) == map_reduced.end())
        map_reduced[cellid] = 0.0;
      if (std::abs(map_reduced[cellid]) < std::abs(val))
        map_reduced[cellid] = val;
    }

  Vector<double> vec(map_reduced.size());
  unsigned int vec_i = 0;
  for (const auto &[cellid, val] : map_reduced) {
    vec[vec_i] = val;
    ++vec_i;
  }
  return vec;
}

double L1NormMap(std::map<CellId, double> &map) {
  Vector<double> vec = GetMPIReducedVectorFromMap(map);
  return vec.l1_norm();
}

double L1NormSizeDividedMap(std::map<CellId, double> &map) {
  Vector<double> vec = GetMPIReducedVectorFromMap(map);
  return vec.l1_norm() / vec.size();
}

double L2NormMap(std::map<CellId, double> &map) {
  Vector<double> vec = GetMPIReducedVectorFromMap(map);
  return vec.l2_norm();
}

double L2NormSizeDividedMap(std::map<CellId, double> &map) {
  Vector<double> vec = GetMPIReducedVectorFromMap(map);
  return vec.l2_norm() / std::sqrt(vec.size());
}

double LInfinityNormMap(std::map<CellId, double> &map) {
  Vector<double> vec = GetMPIReducedVectorFromMap(map);
  return vec.linfty_norm();
}

void L1NormalizeMap(std::map<CellId, double> &map) {
  double l1_norm = L1NormMap(map);
  if (l1_norm == 0)
    throw std::runtime_error("In L1NormalizeMap fucntion, the L2 norm is zero => you have send a "
                             "map with zero values. Check the input map.");

  for (auto &[cellid, val] : map)
    val /= l1_norm;
}

void L1NormalizeSizeIndependentMap(std::map<CellId, double> &map) {
  double l1_norm = L1NormSizeDividedMap(map);
  if (l1_norm == 0)
    throw std::runtime_error(
        "In L1NormalizeSizeIndependentMap fucntion, the L2 norm is zero => you have send a "
        "map with zero values. Check the input map.");

  for (auto &[cellid, val] : map)
    val /= l1_norm;
}

void L2NormalizeMap(std::map<CellId, double> &map) {
  double l2_norm = L2NormMap(map);
  if (l2_norm == 0)
    throw std::runtime_error("In L2NormalizeMap fucntion, the L2 norm is zero => you have send a "
                             "map with zero values. Check the input map.");

  for (auto &[cellid, val] : map)
    val /= l2_norm;
}

void L2NormalizeSizeIndependentMap(std::map<CellId, double> &map) {
  double l2_norm = L2NormSizeDividedMap(map);
  if (l2_norm == 0)
    throw std::runtime_error(
        "In L2NormalizeSizeIndependentMap fucntion, the L2 norm is zero => you have send a "
        "map with zero values. Check the input map.");

  for (auto &[cellid, val] : map)
    val /= l2_norm;
}

void LInfinityNormalizeMap(std::map<CellId, double> &map) {
  double linfinity_norm = LInfinityNormMap(map);
  if (linfinity_norm == 0)
    throw std::runtime_error(
        "In LInfinityNormalizeMap fucntion, the L_infinity norm is zero => you "
        "have send a map with zero values. Check the input map.");

  for (auto &[cellid, val] : map)
    val /= linfinity_norm;
}

double MapDotProduct(const std::map<CellId, double> &map_data1,
                     const std::map<CellId, double> &map_data2) {
  double dot_product = 0;
  for (const auto &[cellid, val] : map_data1)
    dot_product += (val * map_data2.at(cellid));

  dot_product = Utilities::MPI::sum(dot_product, mpi_communicator());

  return dot_product;
}

double SmoothedHeavisideProjection(const double val, const double beta, const double rho_0) {
  return (tanh(beta * rho_0) + tanh(beta * (val - rho_0))) /
         (tanh(beta * rho_0) + tanh(beta * (1 - rho_0)));
}

double SmoothedHeavisideProjectionDer(const double val, const double beta, const double rho_0) {
  return (beta - beta * std::pow(tanh(beta * (val - rho_0)), 2)) /
         (tanh(beta * rho_0) + tanh(beta * (1 - rho_0)));
}

double SetLBeamNotchDensity(const Point<2> &pt, const double notch_radius, const bool antialiasing,
                            const double cell_length) {
  double distance_to_notch_center = pt.distance(Point<2>(0, 0));
  if (antialiasing) {
    if (distance_to_notch_center >= notch_radius + cell_length / 2.0)
      return 1.0; // Fully solid
    else if (distance_to_notch_center <= notch_radius - cell_length / 2.0)
      return 0.0; // Fully void
    else {
      // Linear interpolation for anti-aliasing
      double ratio = (distance_to_notch_center - (notch_radius - cell_length / 2.0)) / cell_length;
      return ratio; // Intermediate value between 0 and 1
    }
  } else {
    return distance_to_notch_center >= notch_radius ? 1.0 : 0.0; // Binary density
  }
}