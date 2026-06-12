#pragma once

#include <LinearElasticMaterial.h>
#include <PNStress.h>

class PNStressLinear : public PNStress {
public:
  PNStressLinear(BVP &bvp, const unsigned int &id);
  ~PNStressLinear() = default;

  // Implement pure virtual methods
  void ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                  double &cell_gradient_density) override;

protected:
  void _ComputeHelpVariables() override;
  void _AssembleAdjointRhs() override;

private:
  std::unique_ptr<MaterialLinElastic> _l_mat;
};