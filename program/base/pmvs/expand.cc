#define _USE_MATH_DEFINES
#include <cmath>

#include "expand.h"
#include "findMatch.h"
#include <algorithm>
#include <iterator>
#include <numeric>

#include "time.h"

using namespace PMVS3;
using namespace std;
using namespace ptch;

Expand::Expand(FindMatch &findMatch) : fm(findMatch) {}

void Expand::Init(void) {}

void Expand::Run(void) {
  fm.count = 0;
  fm.jobs.clear();
  eCounts.resize(fm.CPU);
  fCounts0.resize(fm.CPU);
  fCounts1.resize(fm.CPU);
  pCounts.resize(fm.CPU);
  fill(eCounts.begin(),  eCounts.end(),  0);
  fill(fCounts0.begin(), fCounts0.end(), 0);
  fill(fCounts1.begin(), fCounts1.end(), 0);
  fill(pCounts.begin(),  pCounts.end(),  0);

  time_t starttime = time(NULL);

  fm.po.ClearCounts();
  fm.po.ClearFlags();

  if (!queue.empty()) {
    cerr << "Queue is not empty in expand" << endl;
    exit(1);
  }
  // set queue
  fm.po.CollectPatches(queue);

  cerr << "Expanding patches..." << flush;
  vector<thrd_t> threads(fm.CPU);
  for (int c = 0; c < fm.CPU; ++c)
    thrd_create(&threads[c], &expandThreadTmp, (void *)this);
  for (int c = 0; c < fm.CPU; ++c)
    thrd_join(threads[c], NULL);

  cerr << endl
       << "---- EXPANSION: " << (time(NULL) - starttime) << " secs ----"
       << endl;

  const int trial = accumulate(eCounts.begin(),  eCounts.end(),  0);
  const int fail0 = accumulate(fCounts0.begin(), fCounts0.end(), 0);
  const int fail1 = accumulate(fCounts1.begin(), fCounts1.end(), 0);
  const int pass  = accumulate(pCounts.begin(),  pCounts.end(),  0);
  cerr << "Total pass fail0 fail1 refinepatch: " << trial << ' ' << pass << ' '
       << fail0 << ' ' << fail1 << ' ' << pass + fail1 << endl;
  cerr << "Total pass fail0 fail1 refinepatch: " << 100 * trial / (float)trial
       << ' ' << 100 * pass / (float)trial << ' ' << 100 * fail0 / (float)trial
       << ' ' << 100 * fail1 / (float)trial << ' '
       << 100 * (pass + fail1) / (float)trial << endl;
}

int Expand::expandThreadTmp(void *arg) {
  ((Expand *)arg)->expandThread();
  return 0;
}

void Expand::expandThread(void) {
  mtx_lock(&fm.lock);
  const int id = fm.count++;
  mtx_unlock(&fm.lock);

  while (1) {
    pPatch ppatch;
    int empty = 0;
    mtx_lock(&fm.lock);
    if (queue.empty())
      empty = 1;
    else {
      ppatch = queue.top();
      queue.pop();
    }
    mtx_unlock(&fm.lock);

    if (empty)
      break;

    // For each direction;
    vector<vector<Vec4f>> canCoords;
    findEmptyBlocks(ppatch, canCoords);

    for (int i = 0; i < (int)canCoords.size(); ++i) {
      for (int j = 0; j < (int)canCoords[i].size(); ++j) {
        const int flag = expandSub(ppatch, id, canCoords[i][j]);
        // fail
        if (flag)
          ppatch->dFlag |= (0x0001) << i;
      }
    }
  }
}

void Expand::findEmptyBlocks(const pPatch &ppatch,
                              std::vector<std::vector<Vec4f>> &canCoords) {
  // dnum must be at most 8, because dflag is char
  const int dnum = 6;
  const Patch &patch = *ppatch;

  // Empty six directions
  Vec4f xdir, ydir;
  ortho(ppatch->normal, xdir, ydir);

  // -1: not empty
  // pos: number of free pgrids
  //
  // Check if each direction satisfies both of the following two constraints.
  // a. No neighbor
  // b. At least minImageNumThreshold pgrids without any patches and few
  // counts
  vector<float> fill;
  fill.resize(dnum);
  std::fill(fill.begin(), fill.end(), 0.0f);

  //----------------------------------------------------------------------
  // We look at the effective resolution of each image at the patch.
  // We only use images with good effective resolution to determine
  // empty blocks, because lwo-resolution images can easily satisfy
  // the first condition (neighbors), and no expansion will occur.
  // ----------------------------------------------------------------------
  // Minimum number of images required to obtain high res results, and
  // explor empty blocks.
  const float radius = ComputeRadius(patch);
  const float radiuslow  = radius / 6.0f;  // 2.0f;
  const float radiushigh = radius * 2.5f; // 2.0f;//1.5f;

  vector<pPatch> neighbors;
  fm.po.FindNeighbors(patch, neighbors, 1, 4.0f); // 3.0f);

  vector<pPatch>::iterator bpatch = neighbors.begin();
  vector<pPatch>::iterator epatch = neighbors.end();
  while (bpatch != epatch) {
    const Vec4f diff = (*bpatch)->coord - ppatch->coord;
    Vec2f f2(diff * xdir, diff * ydir);
    const float len = norm(f2);
    if (len < radiuslow || radiushigh < len) {
      ++bpatch;
      continue;
    }

    f2 /= len;
    // unitize(f2);

    float angle = atan2(f2[1], f2[0]);
    if (angle < 0.0)
      angle += 2 * M_PI;

    const float findex = angle / (2 * M_PI / dnum);
    const int   lindex = (int)floor(findex);
    const int   hindex = lindex + 1;

    fill[lindex % dnum] += hindex - findex;
    fill[hindex % dnum] += findex - lindex;
    ++bpatch;
  }

  canCoords.resize(dnum);
  for (int i = 0; i < dnum; ++i) {
    if (0.0f < fill[i])
      // if (0.5f < fill[i])
      continue;

    // If already failed, don't try, because we fail again.
    if (ppatch->dFlag & (0x0001 << i))
      continue;

    const float angle = 2 * M_PI * i / dnum;
    Vec4f canCoord = ppatch->coord + cos(angle) * radius * xdir + sin(angle) * radius * ydir;
    canCoords[i].push_back(canCoord);
  }
}

float Expand::ComputeRadius(const ptch::Patch &patch) {
  const int minnum = 2;
  vector<float> units;
  fm.optim.ComputeUnits(patch, units);
  vector<float> vftmp = units;

  if ((int)units.size() < minnum) {
    cerr << "units size less than minnum: " << (int)units.size() << ' ' << minnum << endl;
    cout << (int)patch.images.size() << endl;
    exit(1);
  }

  nth_element(vftmp.begin(), vftmp.begin() + minnum - 1, vftmp.end());
  // Threshold is the second smallest value with some margin
  // ??? critical
  return (*(vftmp.begin() + minnum - 1)) * fm.cSize;
}

bool Expand::expandSub(const pPatch &orgppatch, const int id, const Vec4f &canCoord) {
  // Choose the closest one
  Patch patch;
  patch.coord  = canCoord;
  patch.normal = orgppatch->normal;
  patch.flag   = 1;

  fm.po.SetGridsImages(patch, orgppatch->images);
  if (patch.images.empty())
    return true;

  //-----------------------------------------------------------------
  // Check bimages and mask. Then, initialize possible visible images
  if (!fm.ps.GetMask(patch.coord, fm.level) || !fm.InsideBimages(patch.coord))
    return true;

  // Check counts and maybe pgrids
  const bool flag = checkCounts(patch);
  if (flag)
    return true;

  // Check edge
  fm.optim.RemoveImagesEdge(patch);
  if (patch.images.empty())
    return true;

  ++eCounts[id];
  //-----------------------------------------------------------------
  // Preprocess
  if (fm.optim.PreProcess(patch, id, 0)) {
    ++fCounts0[id];
    return true;
  }

  //-----------------------------------------------------------------
  fm.optim.RefinePatch(patch, id, 100);

  //-----------------------------------------------------------------
  if (fm.optim.PostProcess(patch, id, 0)) {
    ++fCounts1[id];
    return true;
  }
  ++pCounts[id];

  //-----------------------------------------------------------------
  // Finally
  pPatch ppatch(new Patch(patch));

  // patch.images = orgppatch->images;
  const int add = updateCounts(patch);

  fm.po.AddPatch(ppatch);

  if (add) {
    mtx_lock(&fm.lock);
    queue.push(ppatch);
    mtx_unlock(&fm.lock);
  }

  return false;
}

bool Expand::checkCounts(ptch::Patch &patch) {
  int full  = 0;
  int empty = 0;

  vector<int>::iterator   begin  = patch.images.begin();
  vector<int>::iterator   end    = patch.images.end();
  vector<Vec2i>::iterator begin2 = patch.grids.begin();

  while (begin != end) {
    const int index = *begin;
    if (fm.tNum <= index) {
      ++begin;
      ++begin2;
      continue;
    }

    const int ix = (*begin2)[0];
    const int iy = (*begin2)[1];
    if (ix < 0 || fm.po.gWidths[index] <= ix || iy < 0 || fm.po.gHeights[index] <= iy) {
      ++begin;
      ++begin2;
      continue;
    }

    const int index2 = iy * fm.po.gWidths[index] + ix;

    int flag = 0;
    fm.imageLocks[index].rdlock();
    if (!fm.po.pGrids[index][index2].empty())
      flag = 1;
    fm.imageLocks[index].unlock();
    if (flag) {
      ++full;
      ++begin;
      ++begin2;
      continue;
    }

    // mtx_lock(&fm.countLocks[index]);
    fm.countLocks[index].rdlock();
    if (fm.countThreshold1 <= fm.po.counts[index][index2])
      ++full;
    else
      ++empty;
    //++fm.pos.counts[index][index2];
    fm.countLocks[index].unlock();
    ++begin;
    ++begin2;
  }

  // First expansion is expensive and make the condition strict
  if (fm.depth <= 1) {
    return empty < fm.minImageNumThreshold     && full != 0;
  } else {
    return empty < fm.minImageNumThreshold - 1 && full != 0;
  }
}

bool Expand::updateCounts(const Patch &patch) {
  // Use images and vimages. Loosen when to set add = 1
  int full  = 0;
  int empty = 0;

  {
    vector<int>::const_iterator   begin  = patch.images.begin();
    vector<int>::const_iterator   end    = patch.images.end();
    vector<Vec2i>::const_iterator begin2 = patch.grids.begin();

    while (begin != end) {
      const int index = *begin;
      if (fm.tNum <= index) {
        ++begin;
        ++begin2;
        continue;
      }

      const int ix = (*begin2)[0];
      const int iy = (*begin2)[1];
      if (ix < 0 || fm.po.gWidths[index] <= ix || iy < 0 || fm.po.gHeights[index] <= iy) {
        ++begin;
        ++begin2;
        continue;
      }

      const int index2 = iy * fm.po.gWidths[index] + ix;

      fm.countLocks[index].wrlock();
      if (fm.countThreshold1 <= fm.po.counts[index][index2])
        ++full;
      else
        ++empty;
      ++fm.po.counts[index][index2];

      fm.countLocks[index].unlock();
      ++begin;
      ++begin2;
    }
  }

  {
    vector<int>::const_iterator begin    = patch.vImages.begin();
    vector<int>::const_iterator end      = patch.vImages.end();
    vector<Vec2i>::const_iterator begin2 = patch.vGrids.begin();

    while (begin != end) {
      const int index = *begin;
#ifdef DEBUG
      if (fm.tnum <= index) {
        cerr << "Impossible in updateCounts" << endl;
        exit(1);
      }
#endif

      const int ix = (*begin2)[0];
      const int iy = (*begin2)[1];
      if (ix < 0 || fm.po.gWidths[index] <= ix || iy < 0 || fm.po.gHeights[index] <= iy) {
        ++begin;
        ++begin2;
        continue;
      }

      const int index2 = iy * fm.po.gWidths[index] + ix;

      fm.countLocks[index].wrlock();
      ;
      if (fm.countThreshold1 <= fm.po.counts[index][index2])
        ++full;
      else
        ++empty;
      ++fm.po.counts[index][index2];
      fm.countLocks[index].unlock();
      ++begin;
      ++begin2;
    }
  }

  return empty;
}
