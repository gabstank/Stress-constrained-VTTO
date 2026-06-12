#pragma once

#include <LinearElasticMaterial.h>
#include <ResponseBase.h>

class ComplianceLinear : public ResponseBase {
public:
  ComplianceLinear(BVP &, const unsigned int &id_);

  ~ComplianceLinear();

  double GetFunction() override;
  void ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                  double &cell_gradient_density);

protected:
  std::unique_ptr<MaterialLinElastic> _l_mat;
};