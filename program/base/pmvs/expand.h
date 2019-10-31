#ifndef PMVS3_EXPAND_H
#define PMVS3_EXPAND_H

#include "patchOrganizer.h"
#include <list>
#include <queue>
#include <vector>

namespace PMVS3 {
class FindMatch;

class Expand {
public:
  Expand(FindMatch&findMatch);
  ~Expand(){};

  void Init(void);
  void Run(void);

  float ComputeRadius(const ptch::Patch &patch);

protected:
  bool expandSub(const ptch::pPatch &orgppatch, const int id, const Vec4f &canCoord);

  bool updateCounts(const ptch::Patch &patch);

  bool checkCounts(ptch::Patch &patch);

  void findEmptyBlocks(const ptch::pPatch &ppatch, std::vector<std::vector<Vec4f>> &canCoords);

protected:
  std::priority_queue<ptch::pPatch, std::vector<ptch::pPatch>, P_compare> queue;

  FindMatch &fm;

  //-----------------------------------------------------------------
  // thread related
  //-----------------------------------------------------------------
  void expandThread(void);
  static int expandThreadTmp(void *arg);

  // Number of trials
  std::vector<int> eCounts;
  // Number of failures in the prep
  std::vector<int> fCounts0;
  // Number of failures in the post processing
  std::vector<int> fCounts1;
  // Number passes
  std::vector<int> pCounts;
};
}; // namespace PMVS3

#endif // PMVS3_EXPAND_H
