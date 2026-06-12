// logger.cpp
#include "Logger.h" // Or your actual path
#include <vector> // Ensure included if ProjectParametersType::OptimizationParameters::constraints is a vector

void Logger::writeConvergenceHistory() {

  if (convergence_history_.empty()) {
    std::cerr << "Logger: No data to write." << std::endl;
    return;
  }
  std::string destination_path = parameters().file.destination_path;
  std::string file_name = parameters().file.name;
  std::string base_filename = destination_path + file_name;
  size_t size = 0;
  auto it_iter = convergence_history_.find("iteration");
  if (it_iter != convergence_history_.end() && !it_iter->second.empty()) {
    size = it_iter->second.size();
  } else if (!convergence_history_.empty() && !convergence_history_.begin()->second.empty()) {
    // Fallback if "iteration" is not present or empty, use the size of the first entry
    // This part of logic might need refinement based on guaranteed keys
    size = convergence_history_.begin()->second.size();
    if (size == 0 && it_iter != convergence_history_.end()) {
      // "iteration" key exists but is empty
    } else if (size == 0) {
      std::cerr << "Logger: No data rows to write for convergence history." << std::endl;
      return;
    }
  } else {
    std::cerr << "Logger: No data to determine row count for convergence history." << std::endl;
    return;
  }

  for (const auto &entry : convergence_history_) {
    if (entry.second.size() != size) {
      std::cerr << "Logger Error: Mismatch in data sizes for convergence history." << std::endl;
      // ... (error reporting as before) ...
      return;
    }
  }
  // ... (rest of file writing logic) ...
  std::string logging_filename = base_filename + "-convergence.csv";
  std::ofstream convergence_log_stream(logging_filename, std::ios_base::out);

  if (!convergence_log_stream.is_open()) {
    std::cerr << "Logger Error: Could not open file " << logging_filename << std::endl;
    return;
  }

  bool first_key = true;
  for (const auto &entry : convergence_history_) {
    if (!first_key) {
      convergence_log_stream << " , ";
    }
    convergence_log_stream << entry.first;
    first_key = false;
  }
  convergence_log_stream << std::endl;

  for (size_t i = 0; i < size; ++i) {
    first_key = true;
    for (const auto &entry : convergence_history_) {
      if (!first_key) {
        convergence_log_stream << " , ";
      }
      if (i < entry.second.size()) {
        convergence_log_stream << std::fixed << std::setprecision(10) << entry.second[i];
      }
      first_key = false;
    }
    convergence_log_stream << std::endl;
  }
  convergence_log_stream.close();
  if (this_mpi_process() == 0) { // Assuming you add a way to know MPI rank for messages
    std::cout << "Logger: Convergence history written to " << logging_filename << std::endl;
  }
}

std::string Logger::center(const std::string &s, int w) {
  std::stringstream ss, spaces;
  int pad = w - s.size(); // count excess room to pad
  for (int i = 0; i < pad / 2; ++i)
    spaces << " ";
  ss << spaces.str() << s << spaces.str(); // format with padding
  if (pad > 0 && pad % 2 != 0)             // if pad odd #, add 1 more space
    ss << " ";
  return ss.str();
}

void Logger::printConsoleHeader() {

  if (this_mpi_process() == 0) {
    auto params = parameters();
    std::ios_base::fmtflags original_flags = std::cout.flags();
    char original_fill = std::cout.fill();
    std::streamsize original_precision = std::cout.precision();

    std::cout << std::setprecision(3) << std::scientific << std::setfill(' ');

    int long_width = 11;
    int n_long_entries = 3 + 2 * params.optimization.constraints.size() + 1;
    int width = 2 * n_long_entries + long_width * n_long_entries - 1;

    std::cout << "  ";
    for (int i = 0; i < width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << " |"
              << center("Topology optimization using " + params.optimization.optimizer, width)
              << "|" << std::endl; // Adjusted width for border

    std::cout << " |";
    for (int i = 0; i < width; ++i)
      std::cout << "_";
    std::cout << "|" << std::endl;

    std::cout << std::setfill(' ') << " |" << center("TOTAL", long_width);
    std::cout << std::setfill(' ') << " |" << center("OBJ", long_width) << " |"
              << center("OBJ_RAW", long_width);

    for (unsigned int i = 0; i < params.optimization.constraints.size(); ++i) {
      std::cout << std::setfill(' ') << " |" << center("CON_" + std::to_string(i), long_width)
                << " |" << center("CON_" + std::to_string(i) + "_RAW", long_width);
    }

    std::cout << std::setfill(' ') << " |" << center("CHANGE", long_width) << " |" << std::endl;

    std::cout << " |";
    for (int i = 0; i < width; ++i)
      std::cout << "_";
    std::cout << "|" << std::endl;

    std::cout.flags(original_flags);
    std::cout.fill(original_fill);
    std::cout.precision(original_precision);
  }
}

void Logger::printConsoleDataRow(int opt_iteration, const ResponseValuesType &raw_values,
                                 const ResponseValuesType &mod_values, double change) {

  if (this_mpi_process() == 0) {
    auto params = parameters();
    std::ios_base::fmtflags original_flags = std::cout.flags();
    char original_fill = std::cout.fill();
    std::streamsize original_precision = std::cout.precision();

    std::cout << std::setprecision(3) << std::scientific << std::setfill(' ');
    int long_width = 11;

    // Objective ID - assuming it's fixed or comes from params
    unsigned int obj_id = params.optimization.objective.id; // Use the actual objective ID

    std::cout << " |" << std::setw(long_width) << opt_iteration;

    auto it_mod_obj = mod_values.find(obj_id);
    auto it_raw_obj = raw_values.find(obj_id);

    std::cout << " |" << std::setw(long_width)
              << (it_mod_obj != mod_values.end() ? it_mod_obj->second : 0.0);
    std::cout << " |" << std::setw(long_width)
              << (it_raw_obj != raw_values.end() ? it_raw_obj->second : 0.0);

    for (const auto &constraint_param : params.optimization.constraints) {
      std::stringstream constr_stream, constr_stream_raw;
      auto it_mod_constr = mod_values.find(constraint_param.id);
      auto it_raw_constr = raw_values.find(constraint_param.id);

      constr_stream << std::setprecision(3) << std::scientific
                    << (it_mod_constr != mod_values.end() ? it_mod_constr->second : 0.0);
      constr_stream_raw << std::setprecision(3) << std::scientific
                        << (it_raw_constr != raw_values.end() ? it_raw_constr->second : 0.0);

      std::cout << " |" << std::setw(long_width) << constr_stream.str();
      std::cout << " |" << std::setw(long_width) << constr_stream_raw.str();
    }

    std::cout << " |" << std::setw(long_width) << change;
    std::cout << " |" << std::endl;

    std::cout.flags(original_flags);
    std::cout.fill(original_fill);
    std::cout.precision(original_precision);
  }
}

void Logger::printConsoleFooter(double change) {
  if (this_mpi_process() == 0) {
    auto params = parameters();
    std::cout << std::setprecision(3) << std::scientific << std::setfill(' ');

    int long_width = 11;
    int n_long_entries = 3 + 2 * params.optimization.constraints.size() + 1;
    int width = 2 * n_long_entries + long_width * n_long_entries - 1;

    std::cout << " |";
    for (int i = 0; i < width; ++i)
      std::cout << "_";
    std::cout << "|" << std::endl;

    if (change < params.optimization.convergence.design_change_tol) {
      std::cout << " |" << center("Optimization converged based on design change", width) << "|"
                << std::endl;
    } else {
      std::cout << " |" << center("Optimization iteration limit reached", width) << "|"
                << std::endl;
    }

    std::cout << " |";
    for (int i = 0; i < width; ++i)
      std::cout << "_";
    std::cout << "|" << std::endl;

    std::cout << std::defaultfloat << std::endl;
  }
}