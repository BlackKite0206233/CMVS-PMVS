#ifndef PMVS3_FILTER_H
#define PMVS3_FILTER_H

#include "../numeric/vec2.h"
#include "patch.h"
#include <list>

namespace PMVS3 {

class FindMatch;

class Filter {
public:
  Filter(FindMatch &findMatch);

  void Init(void);
  void Run(void);

  float ComputeGain(const ptch::Patch &patch, const int lock);

  bool FilterQuad(const ptch::Patch &patch, const std::vector<ptch::pPatch> &neighbors) const;

protected:
  void filterOutside(void);
  void filterOutsideThread(void);
  static int filterOutsideThreadTmp(void *arg);

  void filterExact(void);
  void filterExactThread(void);
  static int filterExactThreadTmp(void *arg);

  void filterNeighbor(const int time);
  void filterSmallGroups(void);
  void filterSmallGroupsSub(const int pid, const int id, std::vector<int> &label, std::list<int> &ltmp) const;
  void setDepthMaps(void);
  void setDepthMapsVGridsVPGridsAddPatchV(const int additive);

  void setConf(const int image);

  std::vector<float> gains;

  std::vector<std::vector<int>> newImages, removeImages;
  std::vector<std::vector<TVec2<int>>> newGrids, removeGrids;

  int t_time;
  std::vector<int> rejects;

  //----------------------------------------------------------------------
  // Thread related
  //----------------------------------------------------------------------
  void setDepthMapsThread(void);
  static int setDepthMapsThreadTmp(void *arg);

  void addPatchVThread(void);
  static int addPatchVThreadTmp(void *arg);

  void setVGridsVPGridsThread(void);
  static int setVGridsVPGridsThreadTmp(void *arg);

  void filterNeighborThread(void);
  static int filterNeighborThreadTmp(void *arg);

  FindMatch &fm;
};
}; // namespace PMVS3

#endif // PMVS3_FILTER_H
