#define _USE_MATH_DEFINES
#include <cmath>

#include "findMatch.h"
#include "patchOrganizer.h"
#include <string>

using namespace PMVS3;
using namespace ptch;
using namespace std;

pPatch PatchOrganizer::MAXDEPTH(new   Patch());
pPatch PatchOrganizer::BACKGROUND(new Patch());

PatchOrganizer::PatchOrganizer(FindMatch &findMatch) : fm(findMatch) {}

// change the contents of images from images to indexes
void PatchOrganizer::Image2Index(Patch &patch) {
  // first image has to be tarGet image
  vector<int> newimages;
  for (int i = 0; i < (int)patch.images.size(); ++i) {
    const int index = fm.ps.image2index(patch.images[i]);
    if (index != -1)
      newimages.push_back(index);
  }

  patch.images.swap(newimages);

  // make sure that the reference image is the tarGeting image
  int exist = -1;
  for (int j = 0; j < (int)patch.images.size(); ++j) {
    if (patch.images[j] < fm.tNum) {
      exist = j;
      break;
    }
  }
  if (exist == -1)
    patch.images.clear();
  else if (exist != 0)
    swap(patch.images[0], patch.images[exist]);
}

// change the contents of images from indexes to images
void PatchOrganizer::Index2Image(Patch &patch) {
  for (int i = 0; i < (int)patch.images.size(); ++i)
    patch.images[i] = fm.ps.images[patch.images[i]];
  for (int i = 0; i < (int)patch.vImages.size(); ++i)
    patch.vImages[i] = fm.ps.images[patch.vImages[i]];
}

void PatchOrganizer::Init(void) {
  pGrids.clear();
  pGrids.resize(fm.tNum);
  vpGrids.clear();
  vpGrids.resize(fm.tNum);
  dpGrids.clear();
  dpGrids.resize(fm.tNum);
  counts.clear();
  counts.resize(fm.tNum);

  gWidths.clear();
  gWidths.resize(fm.num);
  gHeights.clear();
  gHeights.resize(fm.num);
  for (int index = 0; index < fm.num; ++index) {
    const int gwidth  = (fm.ps.GetWidth(index,  fm.level) + fm.cSize - 1) / fm.cSize;
    const int gheight = (fm.ps.GetHeight(index, fm.level) + fm.cSize - 1) / fm.cSize;
    gWidths[index]  = gwidth;
    gHeights[index] = gheight;

    if (index < fm.tNum) {
      pGrids[index].resize(gwidth  * gheight);
      vpGrids[index].resize(gwidth * gheight);
      dpGrids[index].resize(gwidth * gheight);
      counts[index].resize(gwidth  * gheight);
      fill(dpGrids[index].begin(), dpGrids[index].end(), MAXDEPTH);
    }
  }
}

void PatchOrganizer::WritePatches2(const std::string prefix, bool bExportPLY, bool bExportptch, bool bExportPSet) {
  CollectPatches(1);

  if (bExportPLY) {
    char buffer[1024];
    sprintf(buffer, "%s.ply", prefix.c_str());
    WritePLY(pPatches, buffer);
  }

  if (bExportptch) {
    char buffer[1024];
    sprintf(buffer, "%s.patch", prefix.c_str());
    ofstream ofstr;
    ofstr.open(buffer);
    ofstr << "PATCHES" << endl << (int)pPatches.size() << endl;
    for (int p = 0; p < (int)pPatches.size(); ++p) {
      Patch patch = *pPatches[p];
      Index2Image(patch);
      ofstr << patch << "\n";
    }
    ofstr.close();
  }

  if (bExportPSet) {
    char buffer[1024];
    sprintf(buffer, "%s.pset", prefix.c_str());
    ofstream ofstr;
    ofstr.open(buffer);
    for (int p = 0; p < (int)pPatches.size(); ++p)
      ofstr << pPatches[p]->coord[0]  << ' ' << pPatches[p]->coord[1]  << ' ' << pPatches[p]->coord[2]  << ' '
            << pPatches[p]->normal[0] << ' ' << pPatches[p]->normal[1] << ' ' << pPatches[p]->normal[2] << "\n";
    ofstr.close();
  }
}

void PatchOrganizer::ReadPatches(void) {
  // Read-in existing reconstructed points. set fix to one for non-tarGeting
  // images
  for (int i = 0; i < fm.tNum; ++i) {
    const int image = fm.images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d", fm.prefix.c_str(), image, fm.level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open())
      continue;

    string header;
    int pnum;
    ifstr >> header >> pnum;
    cerr << image << ' ' << pnum << " patches" << endl;
    for (int p = 0; p < pnum; ++p) {
      pPatch ppatch(new Patch());
      ifstr >> *ppatch;
      ppatch->fix = 0;
      ppatch->vImages.clear();

      Image2Index(*ppatch);
      if (ppatch->images.empty())
        continue;

        // vImages must be tarGeting images
#ifdef DEBUG
      for (int j = 0; j < (int)ppatch->vImages.size(); ++j)
        if (fm.tNum <= ppatch->vImages[j]) {
          cerr
              << "Imposible in readptches. vImages must be tarGeting images"
              << endl
              << "for patches stored in tarGeting images, if visdata2 have "
                 "been consistent"
              << endl;
          exit(1);
        }
#endif
      SetGrids(*ppatch);
      AddPatch(ppatch);
    }

    ifstr.close();
  }

  // For patches in non-tarGeting images
  for (int i = fm.tNum; i < fm.num; ++i) {
    const int image = fm.images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d", fm.prefix.c_str(), image, fm.level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open())
      continue;

    string header;
    int pnum;
    ifstr >> header >> pnum;
    cerr << image << ' ' << pnum << " patches" << endl;

    for (int p = 0; p < pnum; ++p) {
      pPatch ppatch(new Patch());
      ifstr >> *ppatch;
      ppatch->fix = 1;
      ppatch->vImages.clear();

      Image2Index(*ppatch);
      if (ppatch->images.empty())
        continue;

      SetGrids(*ppatch);
      AddPatch(ppatch);
    }

    ifstr.close();
  }
}

void PatchOrganizer::CollectPatches(const bool tarGet) {
  pPatches.clear();

  for (int index = 0; index < fm.tNum; ++index) {
    for (int i = 0; i < (int)pGrids[index].size(); ++i) {
      vector<pPatch>::iterator begin = pGrids[index][i].begin();
      while (begin != pGrids[index][i].end()) {
        (*begin)->id = -1;
        begin++;
      }
    }
  }

  int count = 0;
  for (int index = 0; index < fm.tNum; ++index) {
    for (int i = 0; i < (int)pGrids[index].size(); ++i) {
      vector<pPatch>::iterator begin = pGrids[index][i].begin();
      while (begin != pGrids[index][i].end()) {
        if ((*begin)->id == -1) {
          (*begin)->id = count++;

          if (!tarGet || !(*begin)->fix)
            pPatches.push_back(*begin);
        }
        ++begin;
      }
    }
  }
}

void PatchOrganizer::CollectPatches(std::priority_queue<ptch::pPatch, std::vector<ptch::pPatch>, P_compare> &pqpatches) {
  for (int index = 0; index < fm.tNum; ++index) {
    for (int i = 0; i < (int)pGrids[index].size(); ++i) {
      vector<pPatch>::iterator begin = pGrids[index][i].begin();
      while (begin != pGrids[index][i].end()) {
        if (!(*begin)->flag) {
          (*begin)->flag = 1;
          pqpatches.push(*begin);
        }
        ++begin;
      }
    }
  }
}

void PatchOrganizer::CollectPatches(const int index, std::priority_queue<ptch::pPatch, std::vector<ptch::pPatch>, P_compare> &pqpatches) {
  fm.imageLocks[index].wrlock();
  for (int i = 0; i < (int)pGrids[index].size(); ++i) {
    vector<pPatch>::iterator begin = pGrids[index][i].begin();
    vector<pPatch>::iterator end   = pGrids[index][i].end();

    while (begin != end) {
      if ((*begin)->images[0] == index && !(*begin)->flag) {
        (*begin)->flag = 1;
        pqpatches.push(*begin);
      }
      ++begin;
    }
  }
  fm.imageLocks[index].unlock();
}

// Should be used only for writing
void PatchOrganizer::CollectNonFixPatches(const int index, std::vector<ptch::pPatch> &ppatches) {
  fm.imageLocks[index].wrlock();
  for (int i = 0; i < (int)pGrids[index].size(); ++i) {
    vector<pPatch>::iterator begin = pGrids[index][i].begin();
    vector<pPatch>::iterator end   = pGrids[index][i].end();

    while (begin != end) {
      if ((*begin)->images[0] == index && !(*begin)->fix) {
        ppatches.push_back(*begin);
      }
      ++begin;
    }
  }
  fm.imageLocks[index].unlock();
}

void PatchOrganizer::ClearFlags(void) {
  vector<pPatch>::iterator bppatch = pPatches.begin();
  vector<pPatch>::iterator eppatch = pPatches.end();

  while (bppatch != eppatch) {
    (*bppatch)->flag = 0;
    ++bppatch;
  }
}

void PatchOrganizer::ClearCounts(void) {
  for (int index = 0; index < fm.tNum; ++index) {
    vector<unsigned char>::iterator begin = counts[index].begin();
    vector<unsigned char>::iterator end   = counts[index].end();
    while (begin != end) {
      *begin = (unsigned char)0;
      ++begin;
    }
  }
}

void PatchOrganizer::AddPatch(ptch::pPatch &ppatch) {
  // First handle vImages
  vector<int>::iterator   bimage = ppatch->images.begin();
  vector<int>::iterator   eimage = ppatch->images.end();
  vector<Vec2i>::iterator bgrid  = ppatch->grids.begin();
  while (bimage != eimage) {
    const int index = *bimage;
    if (fm.tNum <= index) {
      ++bimage;
      ++bgrid;
      continue;
    }

    const int index2 = (*bgrid)[1] * gWidths[index] + (*bgrid)[0];
    fm.imageLocks[index].wrlock();
    pGrids[index][index2].push_back(ppatch);
    fm.imageLocks[index].unlock();
    ++bimage;
    ++bgrid;
  }

  // If depth, set vImages
  if (!fm.depth)
    return;

  bimage = ppatch->vImages.begin();
  eimage = ppatch->vImages.end();
  bgrid  = ppatch->vGrids.begin();

  while (bimage != eimage) {
    const int index  = *bimage;
    const int index2 = (*bgrid)[1] * gWidths[index] + (*bgrid)[0];
    fm.imageLocks[index].wrlock();
    vpGrids[index][index2].push_back(ppatch);
    fm.imageLocks[index].unlock();
    ++bimage;
    ++bgrid;
  }

  UpdateDepthMaps(ppatch);
}

void PatchOrganizer::UpdateDepthMaps(pPatch &ppatch) {
  for (int image = 0; image < fm.tNum; ++image) {
    const Vec3f icoord = fm.ps.Project(image, ppatch->coord, fm.level);

    const float fx    = icoord[0] / fm.cSize;
    const int   xs[2] = {(int)floor(fx), (int)ceil(fx)};
    const float fy    = icoord[1] / fm.cSize;
    const int   ys[2] = {(int)floor(fy), (int)ceil(fy)};

    const float depth = fm.ps.photos[image].oAxis * ppatch->coord;

    fm.imageLocks[image].wrlock();
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        if (xs[i] < 0 || gWidths[image] <= xs[i] || ys[j] < 0 || gHeights[image] <= ys[j])
          continue;

        const int index = ys[j] * gWidths[image] + xs[i];
        if (dpGrids[image][index] == MAXDEPTH)
          dpGrids[image][index] = ppatch;
        else {
          const float dtmp = fm.ps.photos[image].oAxis * dpGrids[image][index]->coord;

          if (depth < dtmp)
            dpGrids[image][index] = ppatch;
        }
      }
    }
    fm.imageLocks[image].unlock();
  }
}

void PatchOrganizer::SetGridsImages(ptch::Patch &patch, const std::vector<int> &images) const {
  patch.images.clear();
  patch.grids.clear();
  vector<int>::const_iterator bimage = images.begin();
  vector<int>::const_iterator eimage = images.end();
  while (bimage != eimage) {
    const Vec3f icoord = fm.ps.Project(*bimage, patch.coord, fm.level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / fm.cSize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / fm.cSize;
    if (0 <= ix && ix < gWidths[*bimage] && 
        0 <= iy && iy < gHeights[*bimage]) {
      patch.images.push_back(*bimage);
      patch.grids.push_back(Vec2i(ix, iy));
    }
    ++bimage;
  }
}

void PatchOrganizer::SetGrids(pPatch &ppatch) const { SetGrids(*ppatch); }

void PatchOrganizer::SetGrids(Patch &patch) const {
  patch.grids.clear();
  for (int i = 0; i < (int)patch.images.size(); ++i) {
    const int image = patch.images[i];
    Vec3f icoord = fm.ps.Project(image, patch.coord, fm.level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / fm.cSize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / fm.cSize;
    patch.grids.push_back(TVec2<int>(ix, iy));
  }
}

void PatchOrganizer::SetVImagesVGrids(pPatch &ppatch) {
  SetVImagesVGrids(*ppatch);
}

void PatchOrganizer::SetVImagesVGrids(Patch &patch) {
  vector<int> used;
  used.resize(fm.tNum);
  fill(used.begin(), used.end(), 0);

  vector<int>::iterator bimage = patch.images.begin();
  vector<int>::iterator eimage = patch.images.end();
  while (bimage != eimage) {
    if ((*bimage) < fm.tNum)
      used[*(bimage)] = 1;
    ++bimage;
  }

  bimage = patch.vImages.begin();
  eimage = patch.vImages.end();
  while (bimage != eimage)
    used[*(bimage++)] = 1;

  for (int image = 0; image < fm.tNum; ++image) {
    if (used[image])
      continue;

    int ix, iy;
    if (!IsVisible0(patch, image, ix, iy, fm.neighborThreshold, 1)) {
      continue;
    }

    if (!fm.ps.GetEdge(patch.coord, image, fm.level))
      continue;

    patch.vImages.push_back(image);
    patch.vGrids.push_back(TVec2<int>(ix, iy));
  }
}

void PatchOrganizer::RemovePatch(const pPatch &ppatch) {
  for (int i = 0; i < (int)ppatch->images.size(); ++i) {
    const int image = ppatch->images[i];
    if (fm.tNum <= image)
      continue;

    const int &ix = ppatch->grids[i][0];
    const int &iy = ppatch->grids[i][1];
    const int index = iy * gWidths[image] + ix;
    pGrids[image][index].erase(remove(pGrids[image][index].begin(),
                                        pGrids[image][index].end(), ppatch),
                                 pGrids[image][index].end());
  }

  for (int i = 0; i < (int)ppatch->vImages.size(); ++i) {
    const int image = ppatch->vImages[i];
#ifdef DEBUG
    if (fm.tNum <= image) {
      cerr << "Imposible in removeptch. vImages must be tarGetting images"
           << endl;
      exit(1);
    }
#endif

    const int &ix = ppatch->vGrids[i][0];
    const int &iy = ppatch->vGrids[i][1];
    const int index = iy * gWidths[image] + ix;
    vpGrids[image][index].erase(remove(vpGrids[image][index].begin(),
                                         vpGrids[image][index].end(), ppatch),
                                  vpGrids[image][index].end());
  }
}

bool PatchOrganizer::IsVisible0(const Patch &patch, const int image, int &ix, int &iy, const float strict, const bool lock) {
  const Vec3f icoord = fm.ps.Project(image, patch.coord, fm.level);
  ix = ((int)floor(icoord[0] + 0.5f)) / fm.cSize;
  iy = ((int)floor(icoord[1] + 0.5f)) / fm.cSize;

  return IsVisible(patch, image, ix, iy, strict, lock);
}

bool PatchOrganizer::IsVisible(const Patch &patch, const int image, const int &ix, const int &iy, const float strict, const bool lock) {
  const int &gwidth  = gWidths[image];
  const int &gheight = gHeights[image];

  if (ix < 0 || gwidth <= ix || iy < 0 || gheight <= iy)
    return false;

  if (!fm.depth)
    return true;

  int ans = 0;
  pPatch dppatch  = MAXDEPTH;
  const int index = iy * gwidth + ix;

  if (lock)
    fm.imageLocks[image].rdlock();

  if (dpGrids[image][index] == MAXDEPTH)
    ans = 1;
  else
    dppatch = dpGrids[image][index];

  if (lock)
    fm.imageLocks[image].unlock();

  if (ans)
    return true;

  Vec4f ray = patch.coord - fm.ps.photos[image].center;
  unitize(ray);
  const float diff = ray * (patch.coord - dppatch->coord);
  const float factor = min(2.0, 2.0 + ray * patch.normal);

  return diff < fm.optim.GetUnit(image, patch.coord) * fm.cSize * strict * factor;
}

void PatchOrganizer::FindNeighbors(const ptch::Patch &patch, std::vector<ptch::pPatch> &neighbors, 
                                   const bool lock, const float scale, const int margin, const bool skipvis) {
  const float radius = 1.5 * margin * fm.expand.ComputeRadius(patch);

  vector<int>::const_iterator   bimage = patch.images.begin();
  vector<int>::const_iterator   eimage = patch.images.end();
  vector<Vec2i>::const_iterator bgrid  = patch.grids.begin();
#ifdef DEBUG
  if (patch.images.empty()) {
    cerr << "Empty patches in findCloses" << endl;
    exit(1);
  }
#endif
  float unit = 0.0f;
  for (int i = 0; i < (int)patch.images.size(); ++i)
    unit += fm.optim.GetUnit(patch.images[i], patch.coord);
  unit /= (int)patch.images.size();
  unit *= fm.cSize;

  while (bimage != eimage) {
    if (fm.tNum <= *bimage) {
      ++bimage;
      ++bgrid;
      continue;
    }
    const int image = *bimage;
    const int &ix   = (*bgrid)[0];
    const int &iy   = (*bgrid)[1];
    if (lock)
      fm.imageLocks[image].rdlock();
    for (int j = -margin; j <= margin; ++j) {
      const int ytmp = iy + j;
      if (ytmp < 0 || fm.po.gHeights[image] <= ytmp)
        continue;
      for (int i = -margin; i <= margin; ++i) {
        const int xtmp = ix + i;
        if (xtmp < 0 || fm.po.gWidths[image] <= xtmp)
          continue;
        const int index = ytmp * fm.po.gWidths[image] + xtmp;
        vector<pPatch>::const_iterator bpatch = fm.po.pGrids[image][index].begin();
        vector<pPatch>::const_iterator epatch = fm.po.pGrids[image][index].end();
        while (bpatch != epatch) {
          if (fm.IsNeighborRadius(patch, **bpatch, unit, fm.neighborThreshold * scale, radius))
            neighbors.push_back(*bpatch);
          ++bpatch;
        }
        bpatch = fm.po.vpGrids[image][index].begin();
        epatch = fm.po.vpGrids[image][index].end();
        while (bpatch != epatch) {
          if (fm.IsNeighborRadius(patch, **bpatch, unit, fm.neighborThreshold * scale, radius))
            neighbors.push_back(*bpatch);
          ++bpatch;
        }
      }
    }
    if (lock)
      fm.imageLocks[image].unlock();

    ++bimage;
    ++bgrid;
  }

  if (!skipvis) {
    bimage = patch.vImages.begin();
    eimage = patch.vImages.end();
    bgrid  = patch.vGrids.begin();

    while (bimage != eimage) {
      const int image = *bimage;
      const int &ix   = (*bgrid)[0];
      const int &iy   = (*bgrid)[1];
      if (lock)
        fm.imageLocks[image].rdlock();
      for (int j = -margin; j <= margin; ++j) {
        const int ytmp = iy + j;
        if (ytmp < 0 || fm.po.gHeights[image] <= ytmp)
          continue;
        for (int i = -margin; i <= margin; ++i) {
          const int xtmp = ix + i;
          if (xtmp < 0 || fm.po.gWidths[image] <= xtmp)
            continue;
          const int index = ytmp * fm.po.gWidths[image] + xtmp;
          vector<pPatch>::const_iterator bpatch = fm.po.pGrids[image][index].begin();
          vector<pPatch>::const_iterator epatch = fm.po.pGrids[image][index].end();
          while (bpatch != epatch) {
            if (fm.IsNeighborRadius(patch, **bpatch, unit,  fm.neighborThreshold * scale, radius))
              neighbors.push_back(*bpatch);
            ++bpatch;
          }
          bpatch = fm.po.vpGrids[image][index].begin();
          epatch = fm.po.vpGrids[image][index].end();
          while (bpatch != epatch) {
            if (fm.IsNeighborRadius(patch, **bpatch, unit, fm.neighborThreshold * scale, radius))
              neighbors.push_back(*bpatch);
            ++bpatch;
          }
        }
      }
      if (lock)
        fm.imageLocks[image].unlock();

      ++bimage;
      ++bgrid;
    }

  }

  sort(neighbors.begin(), neighbors.end());
  neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

float PatchOrganizer::ComputeUnit(const ptch::Patch &patch) const {
  float unit = 0.0f;
  for (int i = 0; i < (int)patch.images.size(); ++i)
    unit += fm.optim.GetUnit(patch.images[i], patch.coord);
  unit /= (int)patch.images.size();
  unit *= fm.cSize;
  return unit;
}

void PatchOrganizer::SetScales(ptch::Patch &patch) const {
  const float unit  = fm.optim.GetUnit(patch.images[0], patch.coord);
  const float unit2 = 2.0f * unit;
  Vec4f ray = patch.coord - fm.ps.photos[patch.images[0]].center;
  unitize(ray);

  const int inum = min(fm.tau, (int)patch.images.size());

  // First compute, how many pixel difference per unit along vertical
  // for (int i = 1; i < (int)patch.images.size(); ++i) {
  for (int i = 1; i < inum; ++i) {
    Vec3f diff =
        fm.ps.Project(patch.images[i], patch.coord,               fm.level) -
        fm.ps.Project(patch.images[i], patch.coord - ray * unit2, fm.level);
    patch.dScale += norm(diff);
  }

  // set dScale to the vertical distance where average pixel move is half
  // pixel
  // patch.dScale /= (int)patch.images.size() - 1;
  patch.dScale /= inum - 1;
  patch.dScale  = unit2 / patch.dScale;

  patch.aScale  = atan(patch.dScale / (unit * fm.wSize / 2.0f));
}

// write out results
void PatchOrganizer::WritePLY(const std::vector<pPatch> &patches, const std::string filename) {
  ofstream ofstr;
  ofstr.open(filename.c_str());
  ofstr << "ply" << '\n'
        << "format ascii 1.0" << '\n'
        << "element vertex " << (int)patches.size() << '\n'
        << "property float x" << '\n'
        << "property float y" << '\n'
        << "property float z" << '\n'
        << "property float nx" << '\n'
        << "property float ny" << '\n'
        << "property float nz" << '\n'
        << "property uchar diffuse_red" << '\n'
        << "property uchar diffuse_green" << '\n'
        << "property uchar diffuse_blue" << '\n'
        << "end_header" << '\n';

  vector<pPatch>::const_iterator bpatch = patches.begin();
  vector<pPatch>::const_iterator bend   = patches.end();

  while (bpatch != bend) {
    // Get color
    Vec3i color;

    const int mode = 0;
    // 0: color from images
    // 1: fix
    // 2: angle
    if (mode == 0) {
      int denom = 0;
      Vec3f colorf;
      for (int i = 0; i < (int)(*bpatch)->images.size(); ++i) {
        const int image = (*bpatch)->images[i];
        colorf += fm.ps.GetColor((*bpatch)->coord, image, fm.level);
        denom++;
      }
      colorf /= denom;
      color[0] = min(255, (int)floor(colorf[0] + 0.5f));
      color[1] = min(255, (int)floor(colorf[1] + 0.5f));
      color[2] = min(255, (int)floor(colorf[2] + 0.5f));
    } else if (mode == 1) {
      if ((*bpatch)->tmp == 1.0f) {
        color[0] = 255;
        color[1] = 0;
        color[2] = 0;
      } else {
        color[0] = 255;
        color[1] = 255;
        color[2] = 255;
      }
    } else if (mode == 2) {
      float angle = 0.0f;
      vector<int>::iterator bimage = (*bpatch)->images.begin();
      vector<int>::iterator eimage = (*bpatch)->images.end();

      while (bimage != eimage) {
        const int index = *bimage;
        Vec4f ray = fm.ps.photos[index].center - (*bpatch)->coord;
        ray[3] = 0.0f;
        unitize(ray);

        angle += acos(ray * (*bpatch)->normal);
        ++bimage;
      }

      angle = angle / (M_PI / 2.0f);
      float r, g, b;
      img::Image::Gray2RGB(angle, r, g, b);
      color[0] = (int)(r * 255.0f);
      color[1] = (int)(g * 255.0f);
      color[2] = (int)(b * 255.0f);
    }

    ofstr << (*bpatch)->coord[0]  << ' ' << (*bpatch)->coord[1]  << ' ' << (*bpatch)->coord[2]  << ' ' 
          << (*bpatch)->normal[0] << ' ' << (*bpatch)->normal[1] << ' ' << (*bpatch)->normal[2] << ' '
          << color[0]               << ' ' << color[1]               << ' ' << color[2]               << '\n';
    ++bpatch;
  }
  ofstr.close();
}

void PatchOrganizer::WritePLY(const std::vector<pPatch> &patches, const std::string filename, const std::vector<Vec3i> &colors) {
  ofstream ofstr;
  ofstr.open(filename.c_str());
  ofstr << "ply" << '\n'
        << "format ascii 1.0" << '\n'
        << "element vertex " << (int)patches.size() << '\n'
        << "property float x" << '\n'
        << "property float y" << '\n'
        << "property float z" << '\n'
        << "property float nx" << '\n'
        << "property float ny" << '\n'
        << "property float nz" << '\n'
        << "property uchar diffuse_red" << '\n'
        << "property uchar diffuse_green" << '\n'
        << "property uchar diffuse_blue" << '\n'
        << "end_header" << '\n';

  vector<pPatch>::const_iterator bpatch = patches.begin();
  vector<pPatch>::const_iterator bend   = patches.end();
  vector<Vec3i>::const_iterator  colorb = colors.begin();

  while (bpatch != bend) {
    ofstr << (*bpatch)->coord[0]  << ' ' << (*bpatch)->coord[1]  << ' ' << (*bpatch)->coord[2]  << ' ' 
          << (*bpatch)->normal[0] << ' ' << (*bpatch)->normal[1] << ' ' << (*bpatch)->normal[2] << ' '
          << *colorb << '\n';
    ++bpatch;
    ++colorb;
  }
  ofstr.close();
}
