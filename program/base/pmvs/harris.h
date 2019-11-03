#ifndef PMVS3_HARRIS_H
#define PMVS3_HARRIS_H

#include "../numeric/vec3.h"
#include "detector.h"
#include "point.h"
#include <set>
#include <vector>

namespace PMVS3 {

class Harris : public Detector {
public:
  void Run(const std::vector<unsigned char> &image,
           const std::vector<unsigned char> &mask,
           const std::vector<unsigned char> &edge, 
					 const int width, const int height, const int gspeedup, const float sigma,
           std::multiset<Point> &result);

  virtual ~Harris() {}

protected:
  float sigmaD;
  float sigmaI;

  std::vector<float> gaussD;
  std::vector<float> gaussI;

  std::vector<std::vector<Vec3f>> dIdx;
  std::vector<std::vector<Vec3f>> dIdy;

  std::vector<std::vector<float>> dIdxdIdx;
  std::vector<std::vector<float>> dIdydIdy;
  std::vector<std::vector<float>> dIdxdIdy;

  std::vector<std::vector<float>> response;

  void init(const std::vector<unsigned char> &image, const std::vector<unsigned char> &mask, const std::vector<unsigned char> &edge);

  void setDerivatives(void);
  void preprocess(void);
  void preprocess2(void);
  void setResponse(void);
};
}; // namespace PMVS3

#endif // PMVS3_HARRIS_H
