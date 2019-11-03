#define _USE_MATH_DEFINES
#include <cmath>

#include "detectFeatures.h"
#include "findMatch.h"
#include <ctime>
#include <map>
#include <time.h>

using namespace PMVS3;
using namespace std;
using namespace ptch;

FindMatch::FindMatch(void) : po(*this), seed(*this), expand(*this), filter(*this), optim(*this) {
  debug = 0;
}

FindMatch::~FindMatch() {
  mtx_destroy(&lock);

  for (auto& imgLock : imageLocks)
		imgLock.destroy();
  for (auto& cntLock : countLocks)
		cntLock.destroy();
}

void FindMatch::updateThreshold(void) {
  nccThreshold       -= 0.05f;
  nccThresholdBefore -= 0.05f;

  countThreshold1 = 2;
}

void FindMatch::Init(const Option &option) {
  tImages = option.tImages;
  oImages = option.oImages;
  images.clear();
  images.insert(images.end(), tImages.begin(), tImages.end());
  images.insert(images.end(), oImages.begin(), oImages.end());

  tNum = (int)tImages.size();
  num  = (int)images.size();

  prefix               = option.prefix;
  level                = option.level;
  cSize                = option.cSize;
  nccThreshold         = option.threshold;
  wSize                = option.wSize;
  minImageNumThreshold = option.minImageNum;
  CPU                  = option.CPU;
  setEdge              = option.setEdge;
  sequenceThreshold    = option.sequence;

  junit = 100;
  // This initialization does not matter
  visibleThreshold      = 0.0f;
  visibleThresholdLoose = 0.0f;

  // tau = max(option.minImageNum * 2, min(num, 5));
  tau = min(option.minImageNum * 2, num);

  depth = 0;

  // set target images and other images
  bindexes = option.bindexes;
  visData  = option.visData;
  visData2 = option.visData2;

  //----------------------------------------------------------------------
  mtx_init(&lock, mtx_plain | mtx_recursive);
  imageLocks.resize(num);
  countLocks.resize(num);
  for (int image = 0; image < num; ++image) {
    imageLocks[image].init();
    countLocks[image].init();
  }
  // We set level + 3, to use multi-resolutional texture grabbing
  ps.Init(images, prefix, level + 3, wSize, 1);

  if (setEdge != 0.0f)
    ps.SetEdge(setEdge);
  ps.SetDistances();

  // Detect features if not yet done
  DetectFeatures df;
  const int fcSize = 16;
  df.Run(ps, num, fcSize, level, CPU);

  // Initialize each core member. po should be first
  po.Init();
  seed.Init(df.points);
  expand.Init();
  filter.Init();
  optim.Init();
  //----------------------------------------------------------------------
  // Init thresholds
  angleThreshold0 = 60.0f * M_PI / 180.0f;
  angleThreshold1 = 60.0f * M_PI / 180.0f;

  countThreshold0 = 2;
  countThreshold1 = 4;
  countThreshold2 = 2;

  neighborThreshold  = 0.5f;
  neighborThreshold1 = 1.0f;

  neighborThreshold2 = 1.0f;

  maxAngleThreshold = option.maxAngleThreshold;

  nccThresholdBefore = nccThreshold - 0.3f;

  quadThreshold = option.quadThreshold;

  epThreshold = 2.0f;
}

bool FindMatch::InsideBimages(const Vec4f &coord) const {
  for (const auto& index : bindexes) {
    const Vec3f icoord = ps.Project(index, coord, level);
    if (icoord[0] < 0.0 || ps.GetWidth(index, level)  - 1 < icoord[0] ||
        icoord[1] < 0.0 || ps.GetHeight(index, level) - 1 < icoord[1])
      return false;
  }
  return true;
}

bool FindMatch::IsNeighbor(const ptch::Patch &lhs, const ptch::Patch &rhs, const float neighborThreshold) const {
  const float hunit = (optim.GetUnit(lhs.images[0], lhs.coord) +
                       optim.GetUnit(rhs.images[0], rhs.coord)) /
                      2.0 * cSize;
  return IsNeighbor(lhs, rhs, hunit, neighborThreshold);
}

bool FindMatch::IsNeighbor(const ptch::Patch &lhs, const ptch::Patch &rhs, const float hunit, const float neighborThreshold) const {
  if (lhs.normal * rhs.normal < cos(120.0 * M_PI / 180.0))
    return false;

  const Vec4f diff  = rhs.coord  - lhs.coord;
  const float vunit = lhs.dScale + rhs.dScale;

  const float f0 = lhs.normal * diff;
  const float f1 = rhs.normal * diff;
  float ftmp = (fabs(f0) + fabs(f1)) / 2.0;
  ftmp /= vunit;

  // this may loosen the isneighbor testing. need to tighten (decrease)
  // threshold?
  const float hsize = norm(2 * diff - lhs.normal * f0 - rhs.normal * f1) / 2.0 / hunit;
  if (1.0 < hsize)
    ftmp /= min(2.0f, hsize);

  return ftmp < neighborThreshold;
}

bool FindMatch::IsNeighborRadius(const ptch::Patch &lhs, const ptch::Patch &rhs, const float hunit, const float neighborThreshold, const float radius) const {
  if (lhs.normal * rhs.normal < cos(120.0 * M_PI / 180.0))
    return false;

  const Vec4f diff  = rhs.coord  - lhs.coord;
  const float vunit = lhs.dScale + rhs.dScale;

  const float f0 = lhs.normal * diff;
  const float f1 = rhs.normal * diff;
  float ftmp = (fabs(f0) + fabs(f1)) / 2.0;
  ftmp /= vunit;

  // this may loosen the isneighbor testing. need to tighten (decrease)
  // threshold?
  const float hsize = norm(2 * diff - lhs.normal * f0 - rhs.normal * f1) / 2.0 / hunit;

  // radius check
  if (radius / hunit < hsize)
    return false;

  if (1.0 < hsize)
    ftmp /= min(2.0f, hsize);

  return ftmp < neighborThreshold;
}

void FindMatch::Run(void) {
	clock_t begin = clock();

  //----------------------------------------------------------------------
  // Seed generation
  seed.Run();
  seed.Clear();

  ++depth;
  po.CollectPatches();

  //----------------------------------------------------------------------
  // Expansion
  const int TIME = 3;
  for (int t = 0; t < TIME; ++t) {
    expand.Run();
    filter.Run();

    updateThreshold();

    cout << "STATUS: ";
    for (int i = 0; i < (int)optim.status.size(); ++i) {
      cout << optim.status[i] << ' ';
      if (i % 10 == 9)
        cout << endl;
    }
    cout << endl;

    ++depth;
  }

  cerr << "---- Total: " << (double)(clock() - begin) / CLOCKS_PER_SEC << " secs ----" << endl;
}

void FindMatch::Write(const std::string prefix, const bool bExportPLY, const bool bExportPatch, const bool bExportPSet) {
  po.WritePatches2(prefix, bExportPLY, bExportPatch, bExportPSet);
}
