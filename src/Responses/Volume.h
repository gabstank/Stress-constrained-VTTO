#pragma once

#include <ResponseBase.h>

class Volume : public ResponseBase {

public:
  Volume(BVP &, const unsigned int &);

  ~Volume();

  double GetFunction() override;
  double GetFunctionRef() override;

  /**
   * Computation of the density sensitivity for a given cell.
   *
   * @param cell - cell for which the gradient has to be integrated
   * @param cell_gradient_density - value to be computed
   */
  void ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                  double &cell_gradient_density);

protected:
  bool first_call = true;
  double _solid_volume = 0.;
  double _vol_frac_constraint = 0.;
  // void _OutputResults() const;
};