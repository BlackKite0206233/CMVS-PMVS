#ifndef PMVS3_FINDMATCH_H
#define PMVS3_FINDMATCH_H

#include "patch.h"
#include "rwmutex.h"
#include "tinycthread.h"
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <string>
#include <vector>

#include "../image/photoSet.h"
#include "expand.h"
#include "filter.h"
#include "optim.h"
#include "option.h"
#include "patchOrganizer.h"
#include "seed.h"

namespace PMVS3 {

class FindMatch {
public:
  FindMatch(void);
  virtual ~FindMatch();

  void Init(const PMVS3::Soption &option);
  void Run(void);
  void Write(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet);

  bool InsideBimages(const Vec4f &coord) const;

  bool IsNeighborRadius(const ptch::Patch &lhs, const ptch::Patch &rhs, const float hunit, const float neighborThreshold, const float radius) const;

  bool IsNeighbor(const ptch::Patch &lhs, const ptch::Patch &rhs, const float hunit, const float neighborThreshold) const;
  bool IsNeighbor(const ptch::Patch &lhs, const ptch::Patch &rhs, const float neighborThreshold) const;

  //----------------------------------------------------------------------
  // num of target images
  int tNum;
  // num of total images
  int num;
  // target images
  std::vector<int> tImages;
  // other images where patches are not computed
  std::vector<int> oImages;
  // total images
  std::vector<int> images;

  // prefix
  std::string prefix;
  // level
  int level;
  // cellsize
  int cSize;
  // nccThreshold
  float nccThreshold;
  // windows size
  int wSize;
  // mininum image num threshold
  int minImageNumThreshold;
  // use edge detection or not
  float setEdge;
  // bounding images
  std::vector<int> bindexes;
  // visdata from SfM. m_num x m_num matrix
  std::vector<std::vector<int>> visData;
  // an array of relavant images
  std::vector<std::vector<int>> visData2;
  // sequence Threshold
  int sequenceThreshold;
  // CPU
  int CPU;
  // Threshold on filterQuad
  float quadThreshold;

  // Maximum number of images used in the optimization
  int tau;

  // If patches are dense or not, that is, if we use check(patch) after patch
  // optimization
  int depth;

  //----------------------------------------------------------------------
  // Thresholds
  //----------------------------------------------------------------------
  // For first feature matching. Images within this angle are used in
  // matching.
  float angleThreshold0;
  // tigher angle
  float angleThreshold1;

  // Number of success generation from each seed point
  int countThreshold0;
  // Number of counts, expansion can be tried
  int countThreshold1;

  // Number of trials for each cell in seed
  int countThreshold2;

  // Parameter for isNeighbor in findemptyblocks
  float neighborThreshold;
  // Parameter for isNeighbor in filterOutside
  float neighborThreshold1;
  // Parameter for filterNeighbor
  float neighborThreshold2;

  // ncc threshold before optim
  float nccThresholdBefore;
  // Maximum angle of images must be at least as large as this
  float maxAngleThreshold;

  // visibility consistency threshold
  float visibleThreshold;
  float visibleThresholdLoose;

  // Epipolar distance in seed generation
  float epThreshold;

  //----------------------------------------------------------------------
  // For threads related
  //----------------------------------------------------------------------
  // General lock
  mtx_t lock;
  // For each image
  std::vector<RWMutex> imageLocks;
  std::vector<RWMutex> countLocks;
  // count
  int count;
  // jobs
  std::list<int> jobs;
  // job unit
  int junit;

  //----------------------------------------------------------------------
  // Core members
  //----------------------------------------------------------------------
  // Images
  img::PhotoSet ps;
  // ptch organizer
  PatchOrganizer po;

  // Seed generator
  Seed seed;
  // ptch expansion
  Expand expand;

public:
  // ptch filter
  Filter filter;
  // ptch optimizer
  Optim optim;

  int debug;

protected:
  void init(void);
  void initTargets(void);
  void updateThreshold(void);
  void initImages(void);
};
}; // namespace PMVS3

#endif // PMVS3_FINDMATCH_H
