#include "detector.h"
#include "point.h"
#include <algorithm>
#include <set>

using namespace PMVS3;
using namespace std;

void Detector::SetGaussD(const float sigmaD, std::vector<float> &gaussD) {
  //----------------------------------------------------------------------
  const int marginD = (int)ceil(2 * sigmaD);
  const int sizeD   = 2 * marginD + 1;

  gaussD.resize(sizeD);

  //----------------------------------------------------------------------
  // set gaussD
  float denom = 0.0;
  for (int x = 0; x < sizeD; ++x) {
    int xtmp = x - marginD;
    const float dtmp = xtmp * exp(-(xtmp * xtmp) / (2 * sigmaD * sigmaD));
    gaussD[x] = dtmp;
    if (0.0 < dtmp)
      denom += dtmp;
  }

  for (auto& g : gaussD)
    g /= denom;
}

void Detector::SetGaussI(const float sigmaI, std::vector<float> &gaussI) {
  const int marginI = (int)ceil(2 * sigmaI);
  const int sizeI   = 2 * marginI + 1;

  gaussI.resize(sizeI);

  //----------------------------------------------------------------------
  // set gaussI
  float denom = 0.0;
  for (int x = 0; x < sizeI; ++x) {
    int xtmp = x - marginI;
    const float dtmp = exp(-(xtmp * xtmp) / (2 * sigmaI * sigmaI));
    gaussI[x] = dtmp;
    denom += dtmp;
  }
  
	for (auto& g : gaussI)
		g /= denom;
}

float Detector::setThreshold(std::multiset<Point> &grid) {
  float ave  = 0.0;
  float ave2 = 0.0;
  int count  = 0;

	for (auto& g : grid) {
    count++;
    ave  += g.response;
		ave2 += g.response * g.response;
	}

  if (!count)
    count = 1;

  ave  /= count;
  ave2 /= count;
  ave2  = sqrt(max(0.0f, ave2 - ave * ave));

  // cout << ave << ' ' << ave2 << endl;

  return ave + ave2;
}

bool Detector::isCloseBoundary(const int x, const int y, const int margin) const {
  if (mask.empty())
    return false;

  if (x - margin < 0 || width <= x + margin || y - margin < 0 || height <= y + margin)
    return true;

  for (int j = -margin; j <= margin; ++j) {
    const int ytmp = y + j;
    for (int i = -margin; i <= margin; ++i) {
      const int xtmp = x + i;

      if (!mask[ytmp][xtmp])
        return true;
    }
  }
  return false;
}
