#include <Filter.h>

class DGIProjection : public Filter {
public:
  DGIProjection(Mesh &);
  virtual ~DGIProjection() = default;
  virtual void updateContinuationParameters() override;
  virtual void applyFilter(std::map<CellId, double> &pseudo_densities) override;
  virtual void applyChainrule(std::map<CellId, double> &sensitivities) override;
  std::string getFilterName() const override {
    return FilterNames::DGIProjection;
  } // Override to return name
private:
  void _BuildCellIdCenterMap();
  void _FindNeighbors();
  void _MaxDensity(CellId cell_id, std::map<CellId, double> &reduced_filtered_densities);
  void _MinDensity(CellId cell_id, std::map<CellId, double> &reduced_filtered_densities);

  double _radius;
  double _beta;
  const double _beta_max;
  const double _growth_rate;

  double _min, _max, _mean;

  std::map<CellId, Point<2>> _cellId_center;
  std::map<CellId, double> _cellId_weight;
  std::map<CellId, bool> _cellId_is_bnd;

  std::map<CellId, std::map<CellId, double>> _neighbour_list;
  std::map<CellId, std::map<CellId, bool>> _neighbour_list_double_r;
  std::map<CellId, double> _in_out_map;
  unsigned int _call_counter;

  DoFHandler<2> _dof_handler;
  std::map<CellId, double> _diff_map, _max_map, _min_map, _mean_map, _beta_map, _eta_map;
};