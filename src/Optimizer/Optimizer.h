#pragma once

#include <DataTransfer.h>
#include <GlobalContext.h>

class Optimizer {
public:
  virtual ~Optimizer() = default;

  // This is the core common function. It's a pure virtual function,
  // making Optimizer an abstract class that cannot be instantiated directly.
  virtual void
  DesignUpdate(std::map<dealii::CellId, double> &pseudo_densities,
               const std::map<unsigned int, std::map<dealii::CellId, double>> &gradients,
               const std::map<unsigned int, double> &values) = 0;

  // Provide default implementations if they are not strictly required by all derived classes.
  virtual std::map<std::string, std::map<dealii::CellId, double>> GetDataRelevantForTransfer() {
    return {}; // Default: return empty map
  }
  virtual void SetDataAfterTransfer(const std::map<dealii::CellId, double> & /*data*/,
                                    const std::string & /*name*/) {
    // Default: do nothing
  }
};