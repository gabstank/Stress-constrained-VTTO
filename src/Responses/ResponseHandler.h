#pragma once

// C++ headers
#include <iostream>

// Deal.II headers
#include <deal.II/base/timer.h>

// Project headers
#include <BVP.h>
#include <Miscellaneous.h>
#include <Parameter.h>
#include <ResponseBase.h>

using namespace dealii;

class ResponseHandler {
public:
  ResponseHandler(BVP &);

  // Destructor of base class must always be virtual
  virtual ~ResponseHandler() = default;

  // Getter functionality - just returns the data
  std::map<unsigned int, double> GetRawValues() { return _raw_values; }
  std::map<unsigned int, double> GetRawRefValues() { return _ref_values; }
  std::map<unsigned int, double> GetModifiedValues() { return _mod_values; }
  std::map<unsigned int, std::map<CellId, double>> GetRawGradients() { return _raw_gradients; }
  std::map<unsigned int, std::map<CellId, double>> GetNormalizedGradients() {
    return _normalized_gradients;
  }
  std::map<unsigned int, std::map<CellId, double>> GetChainruleGradients() {
    return _chainrule_gradients;
  }

  void SetChainruleGradients(const std::map<unsigned int, std::map<CellId, double>> &cg) {
    _chainrule_gradients.clear();
    _chainrule_gradients = cg;
  }

  unsigned int GetNResponses() { return _responses.size(); }

  void ComputeResponseValuesGradients();

  void GetAllOutputData(DataOutput<DomainParallelTriaType, DoFHandler<2>> &output) {
    for (auto &[id, response] : _responses)
      response->GetOutputData(output);
  }

private:
  // Adds a response to the response vector
  void _AddResponse(const unsigned int &id, const ResponseNames &response_name);

  // Computing functionality - performs the actual computation
  void _ComputeValues();
  void _ComputeRefValues();
  void _ComputeGradients();
  void _NormalizeValuesGradients();

  void _OutputResults();

  // Response values and gradients
  std::map<unsigned int, double> _raw_values;
  std::map<unsigned int, double> _ref_values;
  std::map<unsigned int, double> _mod_values;
  std::map<unsigned int, std::map<CellId, double>> _raw_gradients;
  std::map<unsigned int, std::map<CellId, double>> _normalized_gradients;
  std::map<unsigned int, std::map<CellId, double>> _chainrule_gradients;

  std::map<unsigned int, std::unique_ptr<ResponseBase>> _responses;
  std::map<unsigned int, std::string> _responses_names;
  std::vector<unsigned int> _response_id;

  BVP &_bvp;

  std::map<unsigned int, double> _scaling_factors;
  bool _scaling_factors_set = false;
  bool _reference_computed = false;
};