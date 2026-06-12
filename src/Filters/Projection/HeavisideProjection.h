#include <Filter.h>
#include <HeavisideKernel.h>

class HeavisideProjection : public Filter {
public:
  HeavisideProjection(Mesh &);
  virtual ~HeavisideProjection() = default;
  virtual void updateContinuationParameters() override;
  virtual void applyFilter(std::map<CellId, double> &pseudo_densities) override;
  virtual void applyChainrule(std::map<CellId, double> &sensitivities) override;
  std::string getFilterName() const override {
    return FilterNames::HeavisideProjection;
  } // Override to return name
private:
  std::unique_ptr<HeavisideKernel> _kernel;
  double _beta;
  const double _beta_max;
  const double _growth_rate;
  std::map<CellId, double> _in_out_map;
};