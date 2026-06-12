#pragma once

// C++ headers
#include <iostream>

// Deal.II headers
#include <deal.II/base/timer.h>

// Project headers
#include <BVP.h>
#include <Miscellaneous.h>
#include <Parameter.h>

using namespace dealii;

// Virtual Base class
class ResponseBase {
public:
  ResponseBase(BVP &, const unsigned int &);

  virtual ~ResponseBase() = default;

  inline unsigned int getID() { return _id; };

  virtual double GetFunction() = 0;
  virtual double GetFunctionRef();

  virtual void ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                          double &cell_gradient_density) {
    (void)cell;
    (void)cell_gradient_density;
    throw std::runtime_error("ComputeCellDensityGradient called from prue "
                             "virtual ResponseBase class.");
  }

  virtual void GetOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {
    (void)output;
    // For most classes there is no need to implement this function, as the base class
    // implementation is empty. Only for the responses that have specific data to output this
    // function needs to be implemented.
  }

protected:
  BVP &_bvp;
  const unsigned int _id;

  double _function_value;
  double _ref_value;
};