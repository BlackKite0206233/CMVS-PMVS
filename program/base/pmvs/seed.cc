#include "seed.h"
#include "findMatch.h"
#include <ctime>
#include <numeric>
#include <time.h>

using namespace img;
using namespace PMVS3;
using namespace ptch;
using namespace std;

Seed::Seed(FindMatch &findMatch) : fm(findMatch) {}

void Seed::Init(const std::vector<std::vector<Point>> &points) {
  pPoints.clear();
  pPoints.resize(fm.num);

  for (int index = 0; index < fm.num; ++index) {
    const int gheight = fm.po.gHeights[index];
    const int gwidth  = fm.po.gWidths[index];
    pPoints[index].resize(gwidth * gheight);
  }

  readPoints(points);
}

void Seed::readPoints(const std::vector<std::vector<Point>> &points) {
  for (int index = 0; index < fm.num; ++index) {
    for (const auto& point : points[index]) {
      pPoint ppoint(new Point(point));
      ppoint->iTmp = index;
      const int ix     = ((int)floor(ppoint->iCoord[0] + 0.5f)) / fm.cSize;
      const int iy     = ((int)floor(ppoint->iCoord[1] + 0.5f)) / fm.cSize;
      const int index2 = iy * fm.po.gWidths[index] + ix;
      pPoints[index][index2].push_back(ppoint);
    }
  }
}

void Seed::Run(void) {
  fm.count = 0;
  fm.jobs.clear();
  sCounts.resize(fm.CPU);
  fCounts0.resize(fm.CPU);
  fCounts1.resize(fm.CPU);
  pCounts.resize(fm.CPU);
  fill(sCounts.begin(),  sCounts.end(),  0);
  fill(fCounts0.begin(), fCounts0.end(), 0);
  fill(fCounts1.begin(), fCounts1.end(), 0);
  fill(pCounts.begin(),  pCounts.end(),  0);

  vector<int> vitmp;
  for (int i = 0; i < fm.tNum; ++i)
    vitmp.push_back(i);

  random_shuffle(vitmp.begin(), vitmp.end());
  fm.jobs.insert(fm.jobs.end(), vitmp.begin(), vitmp.end());

  cerr << "adding seeds " << endl;

  fm.po.ClearCounts();

  // If there already exists a patch, don't use
  for (int index = 0; index < (int)fm.tNum; ++index) 
    for (int j = 0; j < (int)fm.po.pGrids[index].size(); ++j) 
      if (!fm.po.pGrids[index][j].empty())
        fm.po.counts[index][j] = fm.countThreshold2;

	clock_t begin = clock();

  vector<thrd_t> threads(fm.CPU);
  for (auto& t : threads)
    thrd_create(&t, &initialMatchThreadTmp, (void *)this);
  for (auto& t : threads)
    thrd_join(t, NULL);
  //----------------------------------------------------------------------
  cerr << "done" << endl;
  cerr << "---- Initial: " << (double)(clock() - begin) / CLOCKS_PER_SEC << " secs ----" << endl;

  const int trial = accumulate(sCounts.begin(),  sCounts.end(),  0);
  const int fail0 = accumulate(fCounts0.begin(), fCounts0.end(), 0);
  const int fail1 = accumulate(fCounts1.begin(), fCounts1.end(), 0);
  const int pass  = accumulate(pCounts.begin(),  pCounts.end(),  0);

  cerr << "Total pass fail0 fail1 refinepatch: "   << trial << ' ' << pass << ' ' << fail0 << ' ' << fail1 << ' ' << pass + fail1 << endl;
  cerr << "Total pass fail0 fail1 refinepatch: "   << 100 * trial / (float)trial
       << ' ' << 100 * pass           / (float)trial 
       << ' ' << 100 * fail0          / (float)trial
       << ' ' << 100 * fail1          / (float)trial 
       << ' ' << 100 * (pass + fail1) / (float)trial << endl;
}

void Seed::initialMatchThread(void) {
  mtx_lock(&fm.lock);
  const int id = fm.count++;
  mtx_unlock(&fm.lock);

  while (1) {
    int index = -1;
    mtx_lock(&fm.lock);
    if (!fm.jobs.empty()) {
      index = fm.jobs.front();
      fm.jobs.pop_front();
    }
    mtx_unlock(&fm.lock);
    if (index == -1)
      break;

    initialMatch(index, id);
  }
}

int Seed::initialMatchThreadTmp(void *arg) {
  ((Seed *)arg)->initialMatchThread();
  return 0;
}

void Seed::Clear(void) { vector<vector<vector<pPoint>>>().swap(pPoints); }

void Seed::initialMatch(const int index, const int id) {
  vector<int> indexes;
  fm.optim.CollectImages(index, indexes);

  if (fm.tau < (int)indexes.size())
    indexes.resize(fm.tau);

  if (indexes.empty())
    return;

  int totalcount = 0;
  //======================================================================
  // for each feature point, starting from the optical center, keep on
  // matching until we find candidateThreshold patches
  const int gheight = fm.po.gHeights[index];
  const int gwidth  = fm.po.gWidths[index];

  int index2 = -1;
  for (int y = 0; y < gheight; ++y) {
    for (int x = 0; x < gwidth; ++x) {
      ++index2;
      if (!canAdd(index, x, y))
        continue;

      for (const auto& p : pPoints[index][index2]) {
        // collect features that satisfies epipolar geometry
        // constraints and sort them according to the differences of
        // distances between two cameras.
        vector<pPoint> vcp;
        collectCandidates(index, indexes, *p, vcp);

        int count = 0;
        Patch bestpatch;
        //======================================================================
        for (const auto& v : vcp) {
          Patch patch;
          patch.coord  = v->coord;
          patch.normal = fm.ps.photos[index].center - patch.coord;

          unitize(patch.normal);
          patch.normal[3] = 0.0;
          patch.flag = 0;

          ++fm.po.counts[index][index2];
          const int ix     = ((int)floor(v->iCoord[0] + 0.5f)) / fm.cSize;
          const int iy     = ((int)floor(v->iCoord[1] + 0.5f)) / fm.cSize;
          const int index3 = iy * fm.po.gWidths[v->iTmp] + ix;
          if (v->iTmp < fm.tNum)
            ++fm.po.counts[v->iTmp][index3];

          if (!initialMatchSub(index, v->iTmp, id, patch)) {
            ++count;
            if (bestpatch.Score(fm.nccThreshold) < patch.Score(fm.nccThreshold))
              bestpatch = patch;
            if (fm.countThreshold0 <= count)
              break;
          }
        }

        if (count) {
          pPatch ppatch(new Patch(bestpatch));
          fm.po.AddPatch(ppatch);
          ++totalcount;
          break;
        }
      }
    }
  }
  cerr << '(' << index << ',' << totalcount << ')' << flush;
}

void Seed::collectCells(const int index0, const int index1, const Point &p0, std::vector<Vec2i> &cells) {
  Vec3 point(p0.iCoord[0], p0.iCoord[1], p0.iCoord[2]);
#ifdef DEBUG
  if (p0.iCoord[2] != 1.0f) {
    cerr << "Impossible in collectCells" << endl;
    exit(1);
  }
#endif

  Mat3 F;
  img::SetF(fm.ps.photos[index0], fm.ps.photos[index1], F, fm.level);
  const int gwidth  = fm.po.gWidths[index1];
  const int gheight = fm.po.gHeights[index1];

  Vec3 line = transpose(F) * point;
  if (line[0] == 0.0 && line[1] == 0.0) {
    cerr << "Point right on top of the epipole?" << index0 << ' ' << index1 << endl;
    return;
  }
  // vertical
  if (fabs(line[0]) > fabs(line[1])) {
    for (int y = 0; y < gheight; ++y) {
      const float fy = (y + 0.5) * fm.cSize - 0.5f;
      float fx = (-line[1] * fy - line[2]) / line[0];
      fx = max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fx));

      const int ix = ((int)floor(fx + 0.5f)) / fm.cSize;
      if (0 <= ix && ix < gwidth)
        cells.push_back(TVec2<int>(ix, y));
      if (0 <= ix - 1 && ix - 1 < gwidth)
        cells.push_back(TVec2<int>(ix - 1, y));
      if (0 <= ix + 1 && ix + 1 < gwidth)
        cells.push_back(TVec2<int>(ix + 1, y));
    }
  } else {
    for (int x = 0; x < gwidth; ++x) {
      const float fx = (x + 0.5) * fm.cSize - 0.5f;
      float fy = (-line[0] * fx - line[2]) / line[1];
      fy = max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fy));

      const int iy = ((int)floor(fy + 0.5f)) / fm.cSize;
      if (0 <= iy && iy < gheight)
        cells.push_back(TVec2<int>(x, iy));
      if (0 <= iy - 1 && iy - 1 < gheight)
        cells.push_back(TVec2<int>(x, iy - 1));
      if (0 <= iy + 1 && iy + 1 < gheight)
        cells.push_back(TVec2<int>(x, iy + 1));
    }
  }
}

// make sorted array of feature points in images, that satisfy the
// epipolar geometry coming from point in image
void Seed::collectCandidates(const int index, const std::vector<int> &indexes, const Point &point, std::vector<pPoint> &vcp) {
  const Vec3 p0(point.iCoord[0], point.iCoord[1], 1.0);
  for (const auto& indexid: indexes) {
    vector<TVec2<int>> cells;
    collectCells(index, indexid, point, cells);
    Mat3 F;
    img::SetF(fm.ps.photos[index], fm.ps.photos[indexid], F, fm.level);

    for (const auto& c : cells) {
      const int x = c[0];
      const int y = c[1];

      if (!canAdd(indexid, x, y))
        continue;

      const int index2 = y * fm.po.gWidths[indexid] + x;
			for (const auto& p : pPoints[indexid][index2]) {
				if (point.type != p->type) 
					continue;

				const Vec3 p1(p->iCoord[0], p->iCoord[1], 1.0);
				if (fm.epThreshold <= img::ComputeEPD(F, p0, p1)) 
					continue;

				vcp.push_back(p);
			}
    }
  }

  // set distances to response
  vector<pPoint> vcptmp;
  for (auto& v : vcp) {
    unproject(index, v->iTmp, point, *v, v->coord);

    if (fm.ps.photos[index].projection[fm.level][2] * v->coord <= 0.0)
      continue;

    if (!fm.ps.GetMask(v->coord, fm.level) || !fm.InsideBimages(v->coord))
      continue;

    //??? from the closest
    v->response = fabs(norm(v->coord - fm.ps.photos[index].center) - 
                       norm(v->coord - fm.ps.photos[v->iTmp].center));
    vcptmp.push_back(v);
  }

  vcptmp.swap(vcp);
  sort(vcp.begin(), vcp.end());
}

bool Seed::canAdd(const int index, const int x, const int y) {
  if (!fm.ps.GetMask(index, fm.cSize * x, fm.cSize * y, fm.level))
    return false;

  const int index2 = y * fm.po.gWidths[index] + x;

  if (fm.tNum <= index)
    return true;

  // Check if pGrids already contains something
  if (!fm.po.pGrids[index][index2].empty())
    return false;

  //??? critical
  if (fm.countThreshold2 <= fm.po.counts[index][index2])
    return false;

  return true;
}

void Seed::unproject(const int index0, const int index1, const Point &p0, const Point &p1, Vec4f &coord) const {
  Mat4 A;
  A[0][0] = fm.ps.photos[index0].projection[fm.level][0][0] - p0.iCoord[0] * fm.ps.photos[index0].projection[fm.level][2][0];
  A[0][1] = fm.ps.photos[index0].projection[fm.level][0][1] - p0.iCoord[0] * fm.ps.photos[index0].projection[fm.level][2][1];
  A[0][2] = fm.ps.photos[index0].projection[fm.level][0][2] - p0.iCoord[0] * fm.ps.photos[index0].projection[fm.level][2][2];
  A[1][0] = fm.ps.photos[index0].projection[fm.level][1][0] - p0.iCoord[1] * fm.ps.photos[index0].projection[fm.level][2][0];
  A[1][1] = fm.ps.photos[index0].projection[fm.level][1][1] - p0.iCoord[1] * fm.ps.photos[index0].projection[fm.level][2][1];
  A[1][2] = fm.ps.photos[index0].projection[fm.level][1][2] - p0.iCoord[1] * fm.ps.photos[index0].projection[fm.level][2][2];
  A[2][0] = fm.ps.photos[index1].projection[fm.level][0][0] - p1.iCoord[0] * fm.ps.photos[index1].projection[fm.level][2][0];
  A[2][1] = fm.ps.photos[index1].projection[fm.level][0][1] - p1.iCoord[0] * fm.ps.photos[index1].projection[fm.level][2][1];
  A[2][2] = fm.ps.photos[index1].projection[fm.level][0][2] - p1.iCoord[0] * fm.ps.photos[index1].projection[fm.level][2][2];
  A[3][0] = fm.ps.photos[index1].projection[fm.level][1][0] - p1.iCoord[1] * fm.ps.photos[index1].projection[fm.level][2][0];
  A[3][1] = fm.ps.photos[index1].projection[fm.level][1][1] - p1.iCoord[1] * fm.ps.photos[index1].projection[fm.level][2][1];
  A[3][2] = fm.ps.photos[index1].projection[fm.level][1][2] - p1.iCoord[1] * fm.ps.photos[index1].projection[fm.level][2][2];

  Vec4 b;
  b[0] = p0.iCoord[0] * fm.ps.photos[index0].projection[fm.level][2][3] - fm.ps.photos[index0].projection[fm.level][0][3];
  b[1] = p0.iCoord[1] * fm.ps.photos[index0].projection[fm.level][2][3] - fm.ps.photos[index0].projection[fm.level][1][3];
  b[2] = p1.iCoord[0] * fm.ps.photos[index1].projection[fm.level][2][3] - fm.ps.photos[index1].projection[fm.level][0][3];
  b[3] = p1.iCoord[1] * fm.ps.photos[index1].projection[fm.level][2][3] - fm.ps.photos[index1].projection[fm.level][1][3];

  Mat4 AT  = transpose(A);
  Mat4 ATA = AT * A;
  Vec4 ATb = AT * b;

  Mat3 ATA3;
  for (int y = 0; y < 3; ++y)
    for (int x = 0; x < 3; ++x)
      ATA3[y][x] = ATA[y][x];

  Vec3 ATb3;
  for (int y = 0; y < 3; ++y)
    ATb3[y] = ATb[y];

  Mat3 iATA3;
  invert(iATA3, ATA3);
  Vec3 ans = iATA3 * ATb3;
  for (int y = 0; y < 3; ++y)
    coord[y] = ans[y];
  coord[3] = 1.0f;
}

// starting with (index, indexs), set visible images by looking at correlation.
bool Seed::initialMatchSub(const int index0, const int index1, const int id, Patch &patch) {
  //----------------------------------------------------------------------
  patch.images.clear();
  patch.images.push_back(index0);
  patch.images.push_back(index1);

  ++sCounts[id];

  //----------------------------------------------------------------------
  // We know that patch.coord is inside bimages and inside mask
  if (fm.optim.PreProcess(patch, id, 1)) {
    ++fCounts0[id];
    return true;
  }

  //----------------------------------------------------------------------
  fm.optim.RefinePatch(patch, id, 100);

  //----------------------------------------------------------------------
  if (fm.optim.PostProcess(patch, id, 1)) {
    ++fCounts1[id];
    return true;
  }

  ++pCounts[id];
  //----------------------------------------------------------------------
  return false;
}
