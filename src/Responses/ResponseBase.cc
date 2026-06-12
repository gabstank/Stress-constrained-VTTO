#include <ResponseBase.h>

ResponseBase::ResponseBase(BVP &bvp_, const unsigned int &id_) : _bvp(bvp_), _id(id_) {}

double ResponseBase::GetFunctionRef() {
  _ref_value = _function_value;
  return _ref_value;
}