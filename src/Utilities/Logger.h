#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <fstream>
#include <iomanip> // For std::setprecision, std::scientific, std::setfill, std::setw
#include <iostream>
#include <map>
#include <sstream> // For std::stringstream
#include <string>
#include <vector>

#include <Parameter.h>

using ResponseValuesType = std::map<unsigned int, double>;

// Forward declaration if needed by other parts of your project, or include necessary headers
// For example, if parameters() is a global function or part of a globally accessible object
// struct YourParametersType; // Placeholder
// YourParametersType parameters(); // Placeholder
class Logger {
public:
  static Logger &getInstance() {
    static Logger instance;
    return instance;
  }

  void addData(const std::string &key, double value) { convergence_history_[key].push_back(value); }

  void addData(const std::string &key, int value) {
    convergence_history_[key].push_back(static_cast<double>(value));
  }

  void addData(const std::string &key, long long value) {
    convergence_history_[key].push_back(static_cast<double>(value));
  }

  void writeConvergenceHistory(); // Implementation assumed from previous step or in .cpp

  void clearHistory() { convergence_history_.clear(); }

  // New methods for console printing
  void printConsoleHeader();

  void printConsoleDataRow(int opt_iteration, const ResponseValuesType &raw_values,
                           const ResponseValuesType &mod_values, double change);

  void printConsoleFooter(double change);

private:
  Logger() = default;
  ~Logger() = default;
  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;

  std::map<std::string, std::vector<double>> convergence_history_;

  static std::string center(const std::string &input, int width);
};

#endif // LOGGER_H