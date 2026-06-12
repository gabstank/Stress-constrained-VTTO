#include <DesignField.h>
#include <Logger.h>
#include <Volume.h>

Volume::Volume(BVP &bvp_, const unsigned int &id_) : ResponseBase(bvp_, id_) {
  auto params = parameters();

  _vol_frac_constraint = 1.;
  for (const auto &constr : params.optimization.constraints) {
    if (constr.id == this->_id) {
      _vol_frac_constraint = constr.constr_value;
      break;
    }
  }
}

Volume::~Volume() = default;

double Volume::GetFunction() {

  double vol_per_proc = 0.0, full_vol_per_proc;
  auto densities = DesignField::getCurrent();

  for (const auto &cell : this->_bvp._dof_handler.active_cell_iterators()) {

    if (!cell->is_locally_owned() || cell->material_id() == e_padding_mat_id)
      continue;

    if (first_call)
      full_vol_per_proc += cell->measure();
    vol_per_proc += cell->measure() * densities[cell->id()];
  }

  this->_function_value = Utilities::MPI::sum(vol_per_proc, mpi_communicator());
  if (first_call)
    _solid_volume = Utilities::MPI::sum(full_vol_per_proc, mpi_communicator());

  first_call = false;

  auto &logger = Logger::getInstance();
  logger.addData("volume", this->_function_value);

  return this->_function_value;
}

double Volume::GetFunctionRef() {
  this->_ref_value = _solid_volume * _vol_frac_constraint;
  return this->_ref_value;
}

void Volume::ComputeCellDensityGradient(const typename DoFHandler<2>::active_cell_iterator &cell,
                                        double &cell_gradient_density) {
  cell_gradient_density = cell->measure();
}