#pragma once

#include "MMA.h"
#include "Optimizer.h"

class OptimizerFactory {
public:
  static std::unique_ptr<Optimizer> createOptimizer() {
    auto optimizer = parameters().optimization.optimizer;
    if (optimizer == "MMA") {
      return std::make_unique<MMA>();
    } else {
      throw std::runtime_error("Unknown optimizer name in factory: " + optimizer);
    }
  };
};