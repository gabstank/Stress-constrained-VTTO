#pragma once

// C++ headers

// Deal.II headers

// Project headers
#include <BVP.h>
#include <LinearElasticMaterial.h>

class ElasticityLin : public BVP {
public:
  ElasticityLin(Mesh &);
  ~ElasticityLin();
  virtual void Run() override;
  virtual Vector<double> GetPostprocessingData(unsigned int &data_flag) override;

protected:
  void _SetupSystem();
  virtual void _SetupConstraints();
  void _AssembleSystem();
  void _AssembleRhs();
  void _SolveSystem();
  void _SetupQPointHistory();
  void _UpdateQPointHistory();
  void _AssembleSmoothingRhs(unsigned int);
  void _UpdatePostprocessingData(unsigned int);
  void _ComputeRelaxedVMStress();
  void _Postprocess();
  virtual void GetAllOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) override;
  // Getters for mesh adaptivity
  virtual Vector<double> GetRefinementData(const std::string &ref_type) override;
  virtual void _SetElementStresses() override;

  MaterialLinElastic _mat;

  // Struct that stores the information at a specific quadrature point
  struct QPointHistory {
    SymmetricTensor<2, 2> sym_grad_U;
    SymmetricTensor<2, 2> cauchy;
    double vm_stress;
  };
  std::vector<QPointHistory> quadrature_point_history;
  //---------------------------------------------------
  /**
   * @brief Struct to hold all the post processing data.
   */
  struct PostprocessingData {
    Vector<double> domain_vm_stress;
  } postprocessing_data;

}; // end of ElasticityLin class