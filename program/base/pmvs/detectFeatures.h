#ifndef PMVS3_DETECTFEATURES_H
#define PMVS3_DETECTFEATURES_H

/*
 * A main class to detect features
 */

#include "../image/photoSet.h"
#include "point.h"
#include "tinycthread.h"
#include <list>
#include <string>

namespace img {
class PhotoSet;
};

namespace PMVS3 {

class DetectFeatures {
public:
  DetectFeatures(void);
  virtual ~DetectFeatures();

  void Run(const img::PhotoSet &pss, const int num, const int csize, const int level, const int CPU = 1);

  std::vector<std::vector<Point>> points;

protected:
  const img::PhotoSet *ps;
  int cSize;
  int level;

  //----------------------------------------------------------------------
  // thread related
  //----------------------------------------------------------------------
  mtx_t rwLock;
  int CPU;

  std::list<int> jobs;

  void runThread(void);
  static int runThreadTmp(void *arg);
};
}; // namespace PMVS3

#endif // PMVS3_DETECTFEATURES_H
