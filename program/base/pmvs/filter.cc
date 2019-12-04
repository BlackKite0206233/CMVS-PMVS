#include "filter.h"
#include "../numeric/mylapack.h"
#include "findMatch.h"
#include "tinycthread.h"
#include <ctime>
#include <numeric>
#include <time.h>

using namespace ptch;
using namespace PMVS3;
using namespace std;

Filter::Filter(FindMatch &findMatch) : fm(findMatch) {}

void Filter::Init(void) {}

void Filter::Run(void) {
  setDepthMapsVGridsVPGridsAddPatchV(0);

  filterOutside();
  setDepthMapsVGridsVPGridsAddPatchV(1);

  filterExact();
  setDepthMapsVGridsVPGridsAddPatchV(1);

  filterNeighbor(1);
  setDepthMapsVGridsVPGridsAddPatchV(1);

  filterSmallGroups();
  setDepthMapsVGridsVPGridsAddPatchV(1);
}

void Filter::filterOutside(void) {
	clock_t begin = clock();

  cerr << "FilterOutside" << endl;
  //??? notice (1) here to avoid removing fix=1
  fm.po.CollectPatches(1);

  const int psize = (int)fm.po.pPatches.size();
  gains.resize(psize);

  cerr << "mainbody: " << flush;

  fm.count = 0;
  vector<thrd_t> threads(fm.CPU);
  for (auto& t : threads)
    thrd_create(&t, &filterOutsideThreadTmp, (void *)this);
  for (auto& t : threads)
    thrd_join(t, NULL);
  cerr << endl;

  // delete patches with poitive gains
  double ave  = 0.0f;
  double ave2 = 0.0f;
  int denom   = 0;
  int count   = 0;

  for (int p = 0; p < psize; ++p) {
    ave  += gains[p];
    ave2 += gains[p] * gains[p];
    ++denom;

    if (gains[p] < 0.0) {
      fm.po.RemovePatch(fm.po.pPatches[p]);
      count++;
    }
  }

  if (!denom)
    denom = 1;
  ave  /= denom;
  ave2 /= denom;
  ave2  = sqrt(max(0.0, ave2 - ave * ave));
  cerr << "Gain (ave/var): " << ave << ' ' << ave2 << endl;

  cerr << (int)fm.po.pPatches.size() << " -> "
       << (int)fm.po.pPatches.size() - count << " ("
       << 100 * ((int)fm.po.pPatches.size() - count) / (float)fm.po.pPatches.size()
       << "%)\t" << (double)(clock() - begin) / CLOCKS_PER_SEC << " secs" << endl;
}

float Filter::ComputeGain(const ptch::Patch &patch, const bool lock) {
  float gain = patch.Score2(fm.nccThreshold);

  const int size = (int)patch.images.size();
  for (int i = 0; i < size; ++i) {
    const int index = patch.images[i];
    if (fm.tNum <= index)
      continue;

    const int ix     = patch.grids[i][0];
    const int iy     = patch.grids[i][1];
    const int index2 = iy * fm.po.gWidths[index] + ix;

    float maxpressure = 0.0f;
    if (lock)
      fm.imageLocks[index].rdlock();

    for (int j = 0; j < (int)fm.po.pGrids[index][index2].size(); ++j) 
      if (!fm.IsNeighbor(patch, *fm.po.pGrids[index][index2][j], fm.neighborThreshold1))
        maxpressure = max(maxpressure, fm.po.pGrids[index][index2][j]->ncc - fm.nccThreshold);

    if (lock)
      fm.imageLocks[index].unlock();

    gain -= maxpressure;
  }

  const int vsize = (int)patch.vImages.size();
  for (int i = 0; i < vsize; ++i) {
    const int index = patch.vImages[i];
    if (fm.tNum <= index)
      continue;

    const float pdepth = fm.ps.ComputeDepth(index, patch.coord);

    const int ix     = patch.vGrids[i][0];
    const int iy     = patch.vGrids[i][1];
    const int index2 = iy * fm.po.gWidths[index] + ix;
    float maxpressure = 0.0f;

    if (lock)
      fm.imageLocks[index].rdlock();

    for (int j = 0; j < (int)fm.po.pGrids[index][index2].size(); ++j) {
      const float bdepth = fm.ps.ComputeDepth(index, fm.po.pGrids[index][index2][j]->coord);
      if (pdepth < bdepth && !fm.IsNeighbor(patch, *fm.po.pGrids[index][index2][j], fm.neighborThreshold1)) 
        maxpressure = max(maxpressure, fm.po.pGrids[index][index2][j]->ncc - fm.nccThreshold);
    }

    if (lock)
      fm.imageLocks[index].unlock();

    gain -= maxpressure;
  }
  return gain;
}

void Filter::filterOutsideThread(void) {
  mtx_lock(&fm.lock);
  const int id = fm.count++;
  mtx_unlock(&fm.lock);

  const int size  = (int)fm.po.pPatches.size();
  const int itmp  = (int)ceil(size / (float)fm.CPU);
  const int begin = id * itmp;
  const int end   = min(size, (id + 1) * itmp);

  for (int p = begin; p < end; ++p) {
    pPatch &ppatch = fm.po.pPatches[p];
    gains[p] = ppatch->Score2(fm.nccThreshold);

    const int size = (int)ppatch->images.size();
    for (int i = 0; i < size; ++i) {
      const int index = ppatch->images[i];
      if (fm.tNum <= index)
        continue;

      const int ix     = ppatch->grids[i][0];
      const int iy     = ppatch->grids[i][1];
      const int index2 = iy * fm.po.gWidths[index] + ix;

      float maxpressure = 0.0f;
      for (int j = 0; j < (int)fm.po.pGrids[index][index2].size(); ++j) 
        if (!fm.IsNeighbor(*ppatch, *fm.po.pGrids[index][index2][j], fm.neighborThreshold1))
          maxpressure = max(maxpressure, fm.po.pGrids[index][index2][j]->ncc - fm.nccThreshold);

      gains[p] -= maxpressure;
    }

    const int vsize = (int)ppatch->vImages.size();
    for (int i = 0; i < vsize; ++i) {
      const int index = ppatch->vImages[i];
      if (fm.tNum <= index)
        continue;

      const float pdepth = fm.ps.ComputeDepth(index, ppatch->coord);

      const int ix     = ppatch->vGrids[i][0];
      const int iy     = ppatch->vGrids[i][1];
      const int index2 = iy * fm.po.gWidths[index] + ix;
      float maxpressure = 0.0f;

      for (int j = 0; j < (int)fm.po.pGrids[index][index2].size(); ++j) {
        const float bdepth = fm.ps.ComputeDepth( index, fm.po.pGrids[index][index2][j]->coord);
        if (pdepth < bdepth && !fm.IsNeighbor(*ppatch, *fm.po.pGrids[index][index2][j], fm.neighborThreshold1)) 
          maxpressure = max(maxpressure, fm.po.pGrids[index][index2][j]->ncc - fm.nccThreshold);
      }
      gains[p] -= maxpressure;
    }
  }
}

int Filter::filterOutsideThreadTmp(void *arg) {
  ((Filter *)arg)->filterOutsideThread();
  return 0;
}

void Filter::filterExact(void) {
	clock_t begin = clock();
  cerr << "Filter Exact: " << flush;

  //??? cannot use (1) because we use patch.id to set newimages,....
  fm.po.CollectPatches();
  const int psize = (int)fm.po.pPatches.size();

  // dis associate images
  newImages.clear();
  newGrids.clear();
  removeImages.clear();
  removeGrids.clear();
  newImages.resize(psize);
  newGrids.resize(psize);
  removeImages.resize(psize);
  removeGrids.resize(psize);

  fm.count = 0;
  vector<thrd_t> threads(fm.CPU);
  for (auto& t : threads)
    thrd_create(&t, &filterExactThreadTmp, (void *)this);
  for (auto& t : threads)
    thrd_join(t, NULL);
  cerr << endl;

  //----------------------------------------------------------------------
  for (int p = 0; p < psize; ++p) {
    if (fm.po.pPatches[p]->fix)
      continue;

    for (int i = 0; i < (int)removeImages[p].size(); ++i) {
      const int index = removeImages[p][i];
      if (fm.tNum <= index) {
        cerr << "MUST NOT COME HERE" << endl;
        exit(1);
      }
      const int ix     = removeGrids[p][i][0];
      const int iy     = removeGrids[p][i][1];
      const int index2 = iy * fm.po.gWidths[index] + ix;

      fm.po.pGrids[index][index2].erase(
          remove(fm.po.pGrids[index][index2].begin(),
                 fm.po.pGrids[index][index2].end(),
                 fm.po.pPatches[p]),
          fm.po.pGrids[index][index2].end());
    }
  }
  fm.debug = 1;

  int count = 0;
  for (int p = 0; p < psize; ++p) {
    if (fm.po.pPatches[p]->fix)
      continue;

    Patch &patch = *fm.po.pPatches[p];

    // This should be images in targetting images. Has to come before the next
    // for-loop.
    patch.tImages = (int)newImages[p].size();
    for (int i = 0; i < (int)patch.images.size(); ++i) {
      const int index = patch.images[i];
      if (fm.tNum <= index) {
        newImages[p].push_back(patch.images[i]);
        newGrids[p].push_back(patch.grids[i]);
      }
    }
    patch.images.swap(newImages[p]);
    patch.grids.swap(newGrids[p]);

    if (fm.minImageNumThreshold <= (int)patch.images.size()) {
      fm.optim.SetRefImage(patch, 0);
      fm.po.SetGrids(patch);
    }

    if ((int)patch.images.size() < fm.minImageNumThreshold) {
      fm.po.RemovePatch(fm.po.pPatches[p]);
      count++;
    }
  }

  cerr << (int)fm.po.pPatches.size() << " -> "
       << (int)fm.po.pPatches.size() - count << " ("
       << 100 * ((int)fm.po.pPatches.size() - count) / (float)fm.po.pPatches.size()
       << "%)\t" << (double)(clock() - begin) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterExactThread(void) {
  const int psize = (int)fm.po.pPatches.size();
  vector<vector<int>> newimages, removeimages;
  vector<vector<TVec2<int>>> newgrids, removegrids;
  newimages.resize(psize);
  removeimages.resize(psize);
  newgrids.resize(psize);
  removegrids.resize(psize);

  while (1) {
    mtx_lock(&fm.lock);
    const int image = fm.count++;
    mtx_unlock(&fm.lock);

    if (fm.tNum <= image)
      break;

    cerr << '*' << flush;

    const int w = fm.po.gWidths[image];
    const int h = fm.po.gHeights[image];
    int index = -1;
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        ++index;
        for (int i = 0; i < (int)fm.po.pGrids[image][index].size(); ++i) {
          const Patch &patch = *fm.po.pGrids[image][index][i];
          if (patch.fix)
            continue;

          int safe = 0;

          if (fm.po.IsVisible(patch, image, x, y, fm.neighborThreshold1, 0))
            safe = 1;
          // use 4 neighbors?
          else if (0 < x     && fm.po.IsVisible(patch, image, x - 1, y,     fm.neighborThreshold1, 0))
            safe = 1;
          else if (x < w - 1 && fm.po.IsVisible(patch, image, x + 1, y,     fm.neighborThreshold1, 0))
            safe = 1;
          else if (0 < y     && fm.po.IsVisible(patch, image, x,     y - 1, fm.neighborThreshold1, 0))
            safe = 1;
          else if (y < h - 1 && fm.po.IsVisible(patch, image, x,     y + 1, fm.neighborThreshold1, 0))
            safe = 1;

          if (safe) {
            newimages[patch.id].push_back(image);
            newgrids[patch.id].push_back(TVec2<int>(x, y));
          } else {
            removeimages[patch.id].push_back(image);
            removegrids[patch.id].push_back(TVec2<int>(x, y));
          }
        }
      }
    }
  }

  mtx_lock(&fm.lock);
  for (int p = 0; p < psize; ++p) {
		newImages[p].insert(newImages[p].end(), newimages[p].begin(), newimages[p].end());
    newGrids[p].insert( newGrids[p].end(),  newgrids[p].begin(),  newgrids[p].end());
    removeImages[p].insert(removeImages[p].end(), removeimages[p].begin(), removeimages[p].end());
    removeGrids[p].insert( removeGrids[p].end(),  removegrids[p].begin(),  removegrids[p].end());
  }
  mtx_unlock(&fm.lock);
}

int Filter::filterExactThreadTmp(void *arg) {
  ((Filter *)arg)->filterExactThread();
  return 0;
}

void Filter::filterNeighborThread(void) {
  const int size = (int)fm.po.pPatches.size();
  while (1) {
    int jtmp = -1;
    mtx_lock(&fm.lock);
    if (!fm.jobs.empty()) {
      jtmp = fm.jobs.front();
      fm.jobs.pop_front();
    }
    mtx_unlock(&fm.lock);

    if (jtmp == -1)
      break;

    const int begin = fm.junit * jtmp;
    const int end   = min(size, fm.junit * (jtmp + 1));

    for (int p = begin; p < end; ++p) {
      pPatch &ppatch = fm.po.pPatches[p];
      if (rejects[p])
        continue;

      vector<pPatch> neighbors;
      // fm.po.findNeighbors(*ppatch, neighbors, 0, 4, 2);
      fm.po.FindNeighbors(*ppatch, neighbors, 0, 4, 2, 1);

      //?? new filter
      if ((int)neighbors.size() < 6)
        // if ((int)neighbors.size() < 8)
        rejects[p] = t_time + 1;
      else 
        // Fit a quadratic surface
        if (FilterQuad(*ppatch, neighbors))
          rejects[p] = t_time + 1;
    }
  }

  /*
  mtx_lock(&fm.lock);
  const int id = fm.count++;
  mtx_unlock(&fm.lock);

  const int size = (int)fm.po.ppatches.size();
  const int itmp = (int)ceil(size / (float)fm.CPU);
  const int begin = id * itmp;
  const int end = min(size, (id + 1) * itmp);

  for (int p = begin; p < end; ++p) {
    pPatch& ppatch = fm.po.ppatches[p];
    if (rejects[p])
      continue;

    vector<pPatch> neighbors;
    fm.po.findNeighbors(*ppatch, neighbors, 0, 4, 2);

    //?? new filter
    if ((int)neighbors.size() < 6)
    //if ((int)neighbors.size() < 8)
      rejects[p] = time + 1;
    else {
      // Fit a quadratic surface
      if (filterQuad(*ppatch, neighbors))
        rejects[p] = time + 1;
    }
  }
  */
}

bool Filter::FilterQuad(const ptch::Patch &patch, const std::vector<pPatch> &neighbors) const {
  vector<vector<float>> A;
  vector<float> b, x;

  Vec4f xdir, ydir;
  ortho(patch.normal, xdir, ydir);

  const int nsize = (int)neighbors.size();

  float h = 0.0f;
  for (const auto& n : neighbors)
    h += norm(n->coord - patch.coord);
  h /= nsize;

  A.resize(nsize);
  b.resize(nsize);

  vector<float> fxs, fys, fzs;
  fxs.resize(nsize);
  fys.resize(nsize);
  fzs.resize(nsize);
  for (int n = 0; n < nsize; ++n) {
    A[n].resize(5);
    Vec4f diff = neighbors[n]->coord - patch.coord;
    fxs[n] = diff * xdir / h;
    fys[n] = diff * ydir / h;
    fzs[n] = diff * patch.normal;

    A[n][0] = fxs[n] * fxs[n];
    A[n][1] = fys[n] * fys[n];
    A[n][2] = fxs[n] * fys[n];
    A[n][3] = fxs[n];
    A[n][4] = fys[n];
    b[n]    = fzs[n];
  }
  x.resize(5);
  Cmylapack::lls(A, b, x);

  // Compute residual divided by dscale
  const int inum = min(fm.tau, (int)patch.images.size());
  float unit = 0.0;
  // for (int i = 0; i < (int)patch.images.size(); ++i)
  for (const auto& image : patch.images)
    unit += fm.optim.GetUnit(image, patch.coord);
  // unit /= (int)patch.images.size();
  unit /= inum;

  float residual = 0.0f;
  for (int n = 0; n < nsize; ++n) {
    const float res = x[0] * (fxs[n] * fxs[n]) + 
                      x[1] * (fys[n] * fys[n]) +
                      x[2] * (fxs[n] * fys[n]) + 
                      x[3] *  fxs[n] +
                      x[4] *  fys[n] -
														  fzs[n];
    // residual += fabs(res) / neighbors[n]->dscale;
    residual += fabs(res) / unit;
  }

  residual /= (nsize - 5);

  return !(residual < fm.quadThreshold);
}

int Filter::filterNeighborThreadTmp(void *arg) {
  ((Filter *)arg)->filterNeighborThread();
  return 0;
}

void Filter::filterNeighbor(const int times) {
	clock_t begin = clock();

  cerr << "FilterNeighbor:\t" << flush;

  //??? notice (1) to avoid removing fix=1
  fm.po.CollectPatches(1);
  if (fm.po.pPatches.empty())
    return;

  rejects.resize((int)fm.po.pPatches.size());
  fill(rejects.begin(), rejects.end(), 0);

  // Lapack is not thread-safe? Sometimes, the code gets stuck here.
  int count = 0;
  for (t_time = 0; t_time < times; ++t_time) {
    fm.count = 0;

    fm.jobs.clear();
    const int jtmp = (int)ceil(fm.po.pPatches.size() / (float)fm.junit);
    for (int j = 0; j < jtmp; ++j)
      fm.jobs.push_back(j);


    vector<thrd_t> threads(fm.CPU);
    for (auto& t : threads)
      thrd_create(&t, &filterNeighborThreadTmp, (void *)this);
    for (auto& t : threads)
      thrd_join(t, NULL);

    auto breject = rejects.begin();
		for (const auto& patch : fm.po.pPatches) {
			if ((*breject) == t_time + 1) {
				count++;
				fm.po.RemovePatch(patch);
			}

			++breject;
		}
  }

  cerr << (int)fm.po.pPatches.size() << " -> "
       << (int)fm.po.pPatches.size() - count << " ("
       << 100 * ((int)fm.po.pPatches.size() - count) / (float)fm.po.pPatches.size()
       << "%)\t" << (double)(clock() - begin) / CLOCKS_PER_SEC << endl;
}

//----------------------------------------------------------------------
// Take out small connected components
//----------------------------------------------------------------------
void Filter::filterSmallGroups(void) {
	clock_t begin = clock();

  cerr << "FilterGroups:\t" << flush;

  fm.po.CollectPatches();
  if (fm.po.pPatches.empty())
    return;

  const int psize = (int)fm.po.pPatches.size();
  vector<int> label;
  label.resize(psize);
  fill(label.begin(), label.end(), -1);

  list<int> untouch;
  auto bpatch = fm.po.pPatches.begin();
  for (int p = 0; p < psize; ++p, ++bpatch) {
    untouch.push_back(p);
    (*bpatch)->flag = p;
  }

  int id = -1;
  while (!untouch.empty()) {
    const int pid = untouch.front();
    untouch.pop_front();

    if (label[pid] != -1)
      continue;

    label[pid] = ++id;
    list<int> ltmp;
    ltmp.push_back(pid);

    while (!ltmp.empty()) {
      const int ptmp = ltmp.front();
      ltmp.pop_front();

      filterSmallGroupsSub(ptmp, id, label, ltmp);
    }
  }
  id++;

  vector<int> size;
  size.resize(id);
	
	for (const auto& b : label) 
		++size[b];

  const int threshold = max(20, psize / 10000);
  cerr << threshold << endl;

	for (auto& b : size) 
		b = !(b < threshold);

  int count = 0;
  bpatch = fm.po.pPatches.begin();
	for (const auto& b : label) {
		if ((*bpatch)->fix) {
			++bpatch;
			continue;
		}

		if (!size[b]) {
			fm.po.RemovePatch(*bpatch);
			count++;
		}
		++bpatch;
	}

  cerr << (int)fm.po.pPatches.size() << " -> "
       << (int)fm.po.pPatches.size() - count << " ("
       << 100 * ((int)fm.po.pPatches.size() - count) / (float)fm.po.pPatches.size()
       << "%)\t" << (double)(clock() - begin) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterSmallGroupsSub(const int pid, const int id, std::vector<int> &label, std::list<int> &ltmp) const {
  // find neighbors of ptmp and set their ids
  const Patch &patch = *fm.po.pPatches[pid];

  const int index   = patch.images[0];
  const int ix      = patch.grids[0][0];
  const int iy      = patch.grids[0][1];
  const int gwidth  = fm.po.gWidths[index];
  const int gheight = fm.po.gHeights[index];

  for (int y = -1; y <= 1; ++y) {
    const int iytmp = iy + y;
    if (iytmp < 0 || gheight <= iytmp)
      continue;
    for (int x = -1; x <= 1; ++x) {
      const int ixtmp = ix + x;
      if (ixtmp < 0 || gwidth <= ixtmp)
        continue;

      // if (1 < abs(x) + abs(y))
      // continue;

      const int index2 = iytmp * gwidth + ixtmp;

			for (const auto& grid : fm.po.pGrids[index][index2]) {
				const int itmp = grid->flag;
				if (label[itmp] != -1) 
					continue;

				if (fm.IsNeighbor(patch, *grid, fm.neighborThreshold2)) {
					label[itmp] = id;
					ltmp.push_back(itmp);
				}
			}
			for (const auto& grid : fm.po.vpGrids[index][index2]) {
				const int itmp = grid->flag;
				if (label[itmp] != -1) 
					continue;

				if (fm.IsNeighbor(patch, *grid, fm.neighborThreshold2)) {
					label[itmp] = id;
					ltmp.push_back(itmp);
				}
			}
    }
  }
}

void Filter::setDepthMaps(void) {
  // initialize
  for (auto& grid : fm.po.dpGrids) 
    fill(grid.begin(), grid.end(), fm.po.MAXDEPTH);

  fm.count = 0;
  vector<thrd_t> threads(fm.CPU);
  for (auto& t : threads)
    thrd_create(&t, &setDepthMapsThreadTmp, (void *)this);
  for (auto& t : threads)
    thrd_join(t, NULL);
}

int Filter::setDepthMapsThreadTmp(void *arg) {
  ((Filter *)arg)->setDepthMapsThread();
  return 0;
}

void Filter::setDepthMapsThread(void) {
  while (1) {
    mtx_lock(&fm.lock);
    const int index = fm.count++;
    mtx_unlock(&fm.lock);

    if (fm.tNum <= index)
      break;

    const int gwidth  = fm.po.gWidths[index];
    const int gheight = fm.po.gHeights[index];

		for (const auto& patch : fm.po.pPatches) {
			const Vec3f icoord = fm.ps.Project(index, patch->coord, fm.level);
			const float fx     = icoord[0] / fm.cSize;
			const int   xs[2]  = { (int)floor(fx), (int)ceil(fx) };
			const float fy     = icoord[1] / fm.cSize;
			const int   ys[2]  = { (int)floor(fy), (int)ceil(fy) };

			const float depth  = fm.ps.photos[index].oAxis * patch->coord;

			for (int j = 0; j < 2; ++j) {
				for (int i = 0; i < 2; ++i) {
					if (xs[i] < 0 || gwidth <= xs[i] || ys[j] < 0 || gheight <= ys[j])
						continue;
					const int index2 = ys[j] * gwidth + xs[i];

					if (fm.po.dpGrids[index][index2] == fm.po.MAXDEPTH)
						fm.po.dpGrids[index][index2] = patch;
					else 
						if (depth < fm.ps.photos[index].oAxis * fm.po.dpGrids[index][index2]->coord)
							fm.po.dpGrids[index][index2] = patch;
				}
			}
		}
  }
}

void Filter::setDepthMapsVGridsVPGridsAddPatchV(const bool additive) {
  fm.po.CollectPatches();
  setDepthMaps();

  // clear vpGrids
	for (auto& grid : fm.po.vpGrids) 
		for (auto& p : grid) 
			p.clear();

  if (!additive) 
    // initialization
		for (auto& patch : fm.po.pPatches) {
			patch->vImages.clear();
			patch->vGrids.clear();
		}

  fm.count = 0;
  vector<thrd_t> threads0(fm.CPU);
  for (auto& t : threads0)
    thrd_create(&t, &setVGridsVPGridsThreadTmp, (void *)this);
  for (auto& t : threads0)
    thrd_join(t, NULL);

  fm.count = 0;
  vector<thrd_t> threads1(fm.CPU);
  for (auto& t : threads1)
    thrd_create(&t, &addPatchVThreadTmp, (void *)this);
  for (auto& t : threads1)
    thrd_join(t, NULL);
}

int Filter::setVGridsVPGridsThreadTmp(void *arg) {
  ((Filter *)arg)->setVGridsVPGridsThread();
  return 0;
}

int Filter::addPatchVThreadTmp(void *arg) {
  ((Filter *)arg)->addPatchVThread();
  return 0;
}

void Filter::setVGridsVPGridsThread(void) {
  const int noj  = 1000;
  const int size = (int)fm.po.pPatches.size();
  const int job  = max(1, size / (noj - 1));

  while (1) {
    mtx_lock(&fm.lock);
    const int id = fm.count++;
    mtx_unlock(&fm.lock);

    const int begin = id * job;
    const int end   = min(size, (id + 1) * job);

    if (size <= begin)
      break;

    // add patches to vpGrids
    for (int p = begin; p < end; ++p) {
      pPatch &ppatch = fm.po.pPatches[p];
      fm.po.SetVImagesVGrids(ppatch);
    }
  }
}

void Filter::addPatchVThread(void) {
  while (1) {
    mtx_lock(&fm.lock);
    const int index = fm.count++;
    mtx_unlock(&fm.lock);

    if (fm.tNum <= index)
      break;

		for (const auto& patch : fm.po.pPatches) {
			auto bgrid = patch->vGrids.begin();
			for (const auto& image : patch->vImages) {
				if (image == index) {
					const int ix     = (*bgrid)[0];
					const int iy     = (*bgrid)[1];
					const int index2 = iy * fm.po.gWidths[index] + ix;
					fm.po.vpGrids[index][index2].push_back(patch);
					break;
				}

				++bgrid;
			}
		}
  }
}
