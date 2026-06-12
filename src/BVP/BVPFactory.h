#include <BVP.h>
#include <ElasticityLin.h>

// Factory class
class BVPFactory {
public:
  static std::unique_ptr<BVP> createBVP(const std::string &bvp_name, Mesh &mesh) {
    if (bvp_name == "elasticityLin") {
      return std::make_unique<ElasticityLin>(mesh);
    } else {
      throw std::invalid_argument("Invalid BVP type: " + bvp_name);
    }
  }
};