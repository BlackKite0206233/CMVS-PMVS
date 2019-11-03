#define _USE_MATH_DEFINES
#include <cmath>

#include "findMatch.h"
#include "optim.h"
#include <algorithm>
#include <cstdio>
#include <numeric>

#include "nlopt.hpp"

using namespace ptch;
using namespace PMVS3;
using namespace std;

Optim *Optim::one = NULL;

Optim::Optim(FindMatch &findMatch) : fm(findMatch) {
  one = this;

  status.resize(35);
  fill(status.begin(), status.end(), 0);
}

void Optim::Init(void) {
  vect0T.resize(fm.CPU);
  centersT.resize(fm.CPU);
  raysT.resize(fm.CPU);
  indexesT.resize(fm.CPU);
  dScalesT.resize(fm.CPU);
  aScalesT.resize(fm.CPU);
  paramsT.resize(fm.CPU);

  texsT.resize(fm.CPU);
  weightsT.resize(fm.CPU);

  for (int c = 0; c < fm.CPU; ++c) {
    texsT[c].resize(fm.num);
    weightsT[c].resize(fm.num);
    for (auto& t : texsT[c])
      t.resize(3 * fm.wSize * fm.wSize);
  }

  setAxesScales();
}

void Optim::setAxesScales(void) {
  xAxes.resize(fm.num);
  yAxes.resize(fm.num);
  zAxes.resize(fm.num);
  for (int index = 0; index < fm.num; ++index) {
    zAxes[index] = Vec3f(fm.ps.photos[index].oAxis[0],
                         fm.ps.photos[index].oAxis[1],
                         fm.ps.photos[index].oAxis[2]);
    xAxes[index] = Vec3f(fm.ps.photos[index].projection[0][0][0],
                         fm.ps.photos[index].projection[0][0][1],
                         fm.ps.photos[index].projection[0][0][2]);
    yAxes[index] = cross(zAxes[index], xAxes[index]);
    unitize(yAxes[index]);
    xAxes[index] = cross(yAxes[index], zAxes[index]);
  }

  ipScales.resize(fm.num);
  for (int index = 0; index < fm.num; ++index) {
    const Vec4f xaxe(xAxes[index][0], xAxes[index][1], xAxes[index][2], 0.0);
    const Vec4f yaxe(yAxes[index][0], yAxes[index][1], yAxes[index][2], 0.0);

    const float fx  = xaxe * fm.ps.photos[index].projection[0][0];
    const float fy  = yaxe * fm.ps.photos[index].projection[0][1];
    ipScales[index] = fx + fy;
  }
}

void Optim::CollectImages(const int index, std::vector<int> &indexes) const {
  // Find images with constraints angleThreshold, visData,
  // sequenceThreshold, targets. Results are sorted by
  // CphotoSet::distances.
  indexes.clear();
  Vec4f ray0 = fm.ps.photos[index].oAxis;
  ray0[3]    = 0.0f;

  vector<Vec2f> candidates;
  // Search for only related images
  for (const auto& indextmp: fm.visData2[index]) {
    // if (fm.tnum <= indextmp)
    // continue;
    if (fm.sequenceThreshold != -1 && fm.sequenceThreshold < abs(index - indextmp))
      continue;

    Vec4f ray1 = fm.ps.photos[indextmp].oAxis;
    ray1[3]    = 0.0f;

    if (ray0 * ray1 < cos(fm.angleThreshold0))
      continue;

    candidates.push_back(Vec2f(fm.ps.distances[index][indextmp], indextmp));
  }

  sort(candidates.begin(), candidates.end(), Svec2cmp<float>());
  for (int i = 0; i < min(fm.tau, (int)candidates.size()); ++i)
    indexes.push_back((int)candidates[i][1]);
}

bool Optim::PreProcess(Patch &patch, const int id, const int seed) {
  AddImages(patch);

  // Here define reference images, and sort images.
  // Something similar to constraintImages is done inside.
  constraintImages(patch, fm.nccThresholdBefore, id);

  // Fix the reference image and sort the other  tau - 1 images.
  sortImages(patch);

  // Pierre Moulon (it avoid crash in some case)
  if ((int)patch.images.size() > 0) 
    // setSscales should be here to avoid noisy output
    fm.po.SetScales(patch);

  // Check minimum number of images
  if ((int)patch.images.size() < fm.minImageNumThreshold)
    return true;

  if (fm.ps.CheckAngles(patch.coord, patch.images, fm.maxAngleThreshold, fm.angleThreshold1, fm.minImageNumThreshold)) {
    patch.images.clear();
    return true;
  }

  return false;
}

void Optim::filterImagesByAngle(Patch &patch) {
  vector<int> newindexes;

  vector<int>::iterator bimage = patch.images.begin();
  vector<int>::iterator eimage = patch.images.end();

  while (bimage != eimage) {
    const int index = *bimage;
    Vec4f ray = fm.ps.photos[index].center - patch.coord;
    unitize(ray);
    if (ray * patch.normal < cos(fm.angleThreshold1)) {
      // if reference image is deleted, over
      if (bimage == patch.images.begin()) {
        patch.images.clear();
        return;
      }
    } else
      newindexes.push_back(index);
    ++bimage;
  }

  patch.images.swap(newindexes);
}

bool Optim::PostProcess(Patch &patch, const int id, const int seed) {
  if ((int)patch.images.size() < fm.minImageNumThreshold)
    return true;

  if (!fm.ps.GetMask(patch.coord, fm.level) || !fm.InsideBimages(patch.coord))
    return true;

  AddImages(patch);

  constraintImages(patch, fm.nccThreshold, id);
  filterImagesByAngle(patch);

  if ((int)patch.images.size() < fm.minImageNumThreshold)
    return true;

  fm.po.SetGrids(patch);

  SetRefImage(patch, id);
  constraintImages(patch, fm.nccThreshold, id);

  if ((int)patch.images.size() < fm.minImageNumThreshold)
    return true;

  fm.po.SetGrids(patch);

  // set tImages
  patch.tImages = 0;
	for (const auto& image : patch.images)
		if (image < fm.tNum)
			++patch.tImages;

  patch.tmp = patch.Score2(fm.nccThreshold);
  // Set vimages vgrids.
  if (fm.depth) {
    fm.po.SetVImagesVGrids(patch);

    if (2 <= fm.depth && Check(patch))
      return true;
  }
  return false;
}

void Optim::constraintImages(Patch &patch, const float nccThreshold, const int id) {
  vector<float> inccs;
  setINCCs(patch, inccs, patch.images, id, 0);

  //----------------------------------------------------------------------
  // Constraint images
  vector<int> newimages;
  newimages.push_back(patch.images[0]);
  for (int i = 1; i < (int)patch.images.size(); ++i)
    if (inccs[i] < 1.0f - nccThreshold)
      newimages.push_back(patch.images[i]);

  patch.images.swap(newimages);
}

void Optim::SetRefImage(Patch &patch, const int id) {
#ifdef DEBUG
  if (patch.images.empty()) {
    cerr << "empty images" << endl;
    exit(1);
  }
#endif
  //----------------------------------------------------------------------
  // Set the reference image
  // Only for target images
  vector<int> indexes;
	for (const auto& image : patch.images) 
		if (image < fm.tNum) 
			indexes.push_back(image);

  // To avoid segmentation error on alley dataset. (this code is necessary
  // because of the use of filterExact)
  if (indexes.empty()) {
    patch.images.clear();
    return;
  }

  vector<vector<float>> inccs;
  setINCCs(patch, inccs, indexes, id, 1);

  int refindex = -1;
  float refncc = INT_MAX / 2;
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const float sum = accumulate(inccs[i].begin(), inccs[i].end(), 0.0f);
    if (sum < refncc) {
      refncc = sum;
      refindex = i;
    }
  }

  const int refIndex = indexes[refindex];
  for (auto& image : patch.images) {
    if (image == refIndex) {
      const int itmp  = patch.images[0];
      patch.images[0] = refIndex;
			image           = itmp;
      break;
    }
  }
}

// When no sampling was done, this is used
void Optim::setRefConstraintImages(Patch &patch, const float nccThreshold, const int id) {
  //----------------------------------------------------------------------
  // Set the reference image
  vector<vector<float>> inccs;
  setINCCs(patch, inccs, patch.images, id, 1);

  int refindex = -1;
  float refncc = INT_MAX / 2;
  for (int i = 0; i < (int)patch.images.size(); ++i) {
    const float sum = accumulate(inccs[i].begin(), inccs[i].end(), 0.0f);
    if (sum < refncc) {
      refncc   = sum;
      refindex = i;
    }
  }

  // refindex = 0;

  const float robustThreshold = RobustINCC(1.0f - nccThreshold);
  vector<int> newimages;
  newimages.push_back(patch.images[refindex]);
  for (int i = 0; i < (int)patch.images.size(); ++i)
    if (i != refindex && inccs[refindex][i] < robustThreshold)
      newimages.push_back(patch.images[i]);
  patch.images.swap(newimages);
}

void Optim::sortImages(Patch &patch) const {
  const int newm = 1;
  if (newm == 1) {
    const float threshold = 1.0f - cos(10.0 * M_PI / 180.0);
    vector<int> indexes, indexes2;
    vector<float> units, units2;
    vector<Vec4f> rays, rays2;

    ComputeUnits(patch, indexes, units, rays);

    patch.images.clear();
    if (indexes.size() < 2)
      return;

    units[0] = 0.0f;

    while (!indexes.empty()) {
      auto& ite = min_element(units.begin(), units.end());
      const int index = ite - units.begin();

      patch.images.push_back(indexes[index]);

      // Remove other images within 5 degrees
      indexes2.clear();
      units2.clear();
      rays2.clear();
      for (int j = 0; j < (int)rays.size(); ++j) {
        if (j == index)
          continue;

        indexes2.push_back(indexes[j]);
        rays2.push_back(rays[j]);
        const float ftmp = min(threshold, max(threshold / 2.0f, 1.0f - rays[index] * rays[j]));

        units2.push_back(units[j] * (threshold / ftmp));
      }
      indexes2.swap(indexes);
      units2.swap(units);
      rays2.swap(rays);
    }
  } else {
    //----------------------------------------------------------------------
    // Sort and grab the best tau images. All the other images don't
    // matter.  First image is the reference and fixed
    const float threshold = cos(5.0 * M_PI / 180.0);
    vector<int> indexes, indexes2;
    vector<float> units, units2;
    vector<Vec4f> rays, rays2;

    ComputeUnits(patch, indexes, units, rays);

    patch.images.clear();
    if (indexes.size() < 2)
      return;

    units[0] = 0.0f;

    while (!indexes.empty()) {
      // for (int i = 0; i < size; ++i) {
      auto& ite = min_element(units.begin(), units.end());
      const int index = ite - units.begin();

      patch.images.push_back(indexes[index]);

      // Remove other images within 5 degrees
      indexes2.clear();
      units2.clear();
      rays2.clear();
      for (int j = 0; j < (int)rays.size(); ++j) {
        if (rays[index] * rays[j] < threshold) {
          indexes2.push_back(indexes[j]);
          units2.push_back(units[j]);
          rays2.push_back(rays[j]);
        }
      }
      indexes2.swap(indexes);
      units2.swap(units);
      rays2.swap(rays);
    }
  }
}

bool Optim::Check(Patch &patch) {
  const float gain = fm.filter.ComputeGain(patch, 1);
  patch.tmp = gain;

  if (gain < 0.0) {
    patch.images.clear();
    return true;
  }

  {
    vector<pPatch> neighbors;
    fm.po.FindNeighbors(patch, neighbors, 1, 4, 2);
    // Only check when enough number of neighbors
    if (6 < (int)neighbors.size() &&
        // if (8 < (int)neighbors.size() &&
        fm.filter.FilterQuad(patch, neighbors)) {
      patch.images.clear();
      return true;
    }
  }

  return false;
}

void Optim::RemoveImagesEdge(ptch::Patch &patch) const {
  vector<int> newindexes;

	for (const auto& image : patch.images)
		if (fm.ps.GetEdge(patch.coord, image, fm.level))
			newindexes.push_back(image);

  patch.images.swap(newindexes);
}

void Optim::AddImages(ptch::Patch &patch) const {
  // take into account edge
  vector<int> used;
  used.resize(fm.num);
  for (auto& u : used)
    u = 0;

	for (const auto& image : patch.images) {
		used[image] = 1;
	}

  const float athreshold = cos(fm.angleThreshold0);
	for (const auto& image : fm.visData2[patch.images[0]]) {
		if (used[image])
			continue;

		const Vec3f icoord = fm.ps.Project(image, patch.coord, fm.level);
		if (icoord[0] < 0.0f || fm.ps.GetWidth(image,  fm.level) - 1 <= icoord[0] ||
			  icoord[1] < 0.0f || fm.ps.GetHeight(image, fm.level) - 1 <= icoord[1]) 
			continue;

		if (!fm.ps.GetEdge(patch.coord, image, fm.level)) 
			continue;

		Vec4f ray = fm.ps.photos[image].center - patch.coord;
		unitize(ray);
		const float ftmp = ray * patch.normal;

		if (athreshold <= ftmp)
			patch.images.push_back(image);
	}
}

void Optim::ComputeUnits(const ptch::Patch &patch, std::vector<float> &units) const {
  const int size = (int)patch.images.size();
  units.resize(size);

  auto& bfine = units.begin();

	for (const auto& image : patch.images) {
		*bfine = INT_MAX / 2;

		*bfine = GetUnit(image, patch.coord);
		Vec4f ray = fm.ps.photos[image].center - patch.coord;
		unitize(ray);
		const float denom = ray * patch.normal;
		if (0.0 < denom)
			*bfine /= denom;

		++bfine;
	}
}

void Optim::ComputeUnits(const ptch::Patch &patch, std::vector<int> &indexes, std::vector<float> &units, std::vector<Vec4f> &rays) const {
	for (const auto& image : patch.images) {
		Vec4f ray = fm.ps.photos[image].center - patch.coord;
		unitize(ray);
		const float dot = ray * patch.normal;
		if (dot <= 0.0f)
			continue;

		const float scale = GetUnit(image, patch.coord);
		const float fine  = scale / dot;

		indexes.push_back(image);
		units.push_back(fine);
		rays.push_back(ray);
	}
}

void Optim::RefinePatch(Patch &patch, const int id, const int time) {
  if (!RefinePatchBFGS(patch, id, 1000, 1))
    std::cout << "refinePatchBFGS failed!" << std::endl;

  if (patch.images.empty())
    return;
}

//----------------------------------------------------------------------
// BFGS functions
//----------------------------------------------------------------------

double Optim::my_f(unsigned n, const double *x, double *grad, void *my_func_data) {
  double xs[3] = {x[0], x[1], x[2]};
  const int id = *((int *)my_func_data);

  const float angle1 = xs[1] * one->aScalesT[id];
  const float angle2 = xs[2] * one->aScalesT[id];

  double ret = 0.0;

  //?????
  const double bias = 0.0f; // 2.0 - exp(- angle1 * angle1 / sigma2) - exp(-
                            // angle2 * angle2 / sigma2);

  Vec4f coord, normal;
  one->decode(coord, normal, xs, id);

  const int index = one->indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  one->GetPAxes(index, coord, normal, pxaxis, pyaxis);

  const int size    = min(one->fm.tau, (int)one->indexesT[id].size());
  const int mininum = min(one->fm.minImageNumThreshold, size);

  for (int i = 0; i < size; ++i) {
    if (!one->grabTex(coord, pxaxis, pyaxis, normal, one->indexesT[id][i], one->fm.wSize, one->texsT[id][i]))
      one->Normalize(one->texsT[id][i]);
  }

  const int pairwise = 0;
  if (pairwise) {
    double ans = 0.0f;
    int denom  = 0;
    for (int i = 0; i < size; ++i) {
      for (int j = i + 1; j < size; ++j) {
        if (one->texsT[id][i].empty() || one->texsT[id][j].empty())
          continue;

        ans += RobustINCC(1.0 - one->Dot(one->texsT[id][i], one->texsT[id][j]));
        denom++;
      }
    }
    if (denom <
        // one->fm.minImageNumThreshold *
        //(one->fm.minImageNumThreshold - 1) / 2)
        mininum * (mininum - 1) / 2)
      ret = 2.0f;
    else
      ret = ans / denom + bias;
  } else {
    if (one->texsT[id][0].empty())
      return 2.0;

    double ans = 0.0f;
    int denom  = 0;
    for (int i = 1; i < size; ++i) {
      if (one->texsT[id][i].empty())
        continue;
      ans += RobustINCC(1.0 - one->Dot(one->texsT[id][0], one->texsT[id][i]));
      denom++;
    }
    // if (denom < one->fm.minImageNumThreshold - 1)
    ret = (denom < mininum - 1) ? 2.0f : ans / denom + bias;
  }

  return ret;
}

bool Optim::RefinePatchBFGS(Patch &patch, const int id, const int time, const int ncc) {
  int idtmp = id;

  centersT[id] = patch.coord;
  raysT[id]    = patch.coord - fm.ps.photos[patch.images[0]].center;
  unitize(raysT[id]);
  indexesT[id] = patch.images;

  dScalesT[id] = patch.dScale;
  aScalesT[id] = M_PI / 48.0f; // patch.ascale;

  ComputeUnits(patch, weightsT[id]);
  for (auto& w : weightsT[id])
    w = min(1.0f, weightsT[id][0] / w);
  weightsT[id][0] = 1.0f;

  double p[3];
  encode(patch.coord, patch.normal, p, id);

  double min_angle = -23.99999; //(- M_PI / 2.0) / one->ascalesT[id];
  double max_angle =  23.99999;  //(M_PI / 2.0) / one->ascalesT[id];

  std::vector<double> lower_bounds(3);
  lower_bounds[0] = -HUGE_VAL; // Not bound
  lower_bounds[1] = min_angle;
  lower_bounds[2] = min_angle;
  std::vector<double> upper_bounds(3);
  upper_bounds[0] = HUGE_VAL; // Not bound
  upper_bounds[1] = max_angle;
  upper_bounds[2] = max_angle;

  bool success = false;

  try {
    // LN_NELDERMEAD: Corresponds to the N-Simplex-Algorithm of GSL, that was
    // used originally here LN_SBPLX LN_COBYLA LN_BOBYQA LN_PRAXIS
    nlopt::opt opt(nlopt::LN_BOBYQA, 3);
    opt.set_min_objective(my_f, &idtmp);
    opt.set_xtol_rel(1.e-7);
    opt.set_maxeval(time);

    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);

    std::vector<double> x(3);
    for (int i = 0; i < 3; i++) 
      // NLOPT returns an error if x is not within the bounds
      x[i] = max(min(p[i], upper_bounds[i]), lower_bounds[i]);

    double minf;
    nlopt::srand(1);
    nlopt::result result = opt.optimize(x, minf);

    p[0] = x[0];
    p[1] = x[1];
    p[2] = x[2];

    success = (result == nlopt::SUCCESS      || result == nlopt::STOPVAL_REACHED ||
               result == nlopt::FTOL_REACHED || result == nlopt::XTOL_REACHED);
  } catch (std::exception &e) {
    success = false;
  }

  if (success) {
    decode(patch.coord, patch.normal, p, id);

    patch.ncc = 1.0 - UnrobustINCC(ComputeINCC(patch.coord, patch.normal, patch.images, id, 1));
  } else 
    return false;

  return true;
}

void Optim::encode(const Vec4f &coord, double *const vect, const int id) const {
  vect[0] = (coord - centersT[id]) * raysT[id] / dScalesT[id];
}

void Optim::encode(const Vec4f &coord, const Vec4f &normal, double *const vect, const int id) const {
  encode(coord, vect, id);

  const int   image = indexesT[id][0];
  const float fx    = xAxes[image] * proj(normal); // projects from 4D to 3D, divide by last value
  const float fy    = yAxes[image] * proj(normal);
  const float fz    = zAxes[image] * proj(normal);

  vect[2] = asin(max(-1.0f, min(1.0f, fy)));
  const float cosb = cos(vect[2]);

  if (cosb == 0.0)
    vect[1] = 0.0;
  else {
    const float sina = fx / cosb;
    const float cosa = -fz / cosb;
    vect[1] = acos(max(-1.0f, min(1.0f, cosa)));
    if (sina < 0.0)
      vect[1] = -vect[1];
  }

  vect[1] = vect[1] / aScalesT[id];
  vect[2] = vect[2] / aScalesT[id];
}

void Optim::decode(Vec4f &coord, Vec4f &normal, const double *const vect, const int id) const {
  decode(coord, vect, id);
  const int image = indexesT[id][0];

  const float angle1 = vect[1] * aScalesT[id];
  const float angle2 = vect[2] * aScalesT[id];

  const float fx =  sin(angle1) * cos(angle2);
  const float fy =  sin(angle2);
  const float fz = -cos(angle1) * cos(angle2);

  Vec3f ftmp = xAxes[image] * fx + yAxes[image] * fy + zAxes[image] * fz;
  normal = Vec4f(ftmp[0], ftmp[1], ftmp[2], 0.0f);
}

void Optim::decode(Vec4f &coord, const double *const vect, const int id) const {
  coord = centersT[id] + dScalesT[id] * vect[0] * raysT[id];
}

void Optim::setINCCs(const ptch::Patch &patch, std::vector<float> &inccs, const std::vector<int> &indexes, const int id, const bool robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  GetPAxes(index, patch.coord, patch.normal, pxaxis, pyaxis);

  vector<vector<float>> &texs = texsT[id];

  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) 
    if (!grabTex(patch.coord, pxaxis, pyaxis, patch.normal, indexes[i], fm.wSize, texs[i])) 
      Normalize(texs[i]);

  inccs.resize(size);
  if (texs[0].empty()) {
    fill(inccs.begin(), inccs.end(), 2.0f);
    return;
  }

  for (int i = 0; i < size; ++i) {
    if (!i)
      inccs[i] = 0.0f;
    else if (!texs[i].empty())
      inccs[i] = (!robust) ? 1.0f - Dot(texs[0], texs[i]) : RobustINCC(1.0f - Dot(texs[0], texs[i]));
    else
      inccs[i] = 2.0f;
  }
}

void Optim::setINCCs(const ptch::Patch &patch, std::vector<std::vector<float>> &inccs, const std::vector<int> &indexes, const int id, const bool robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  GetPAxes(index, patch.coord, patch.normal, pxaxis, pyaxis);

  vector<vector<float>> &texs = texsT[id];

  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) 
    if (!grabTex(patch.coord, pxaxis, pyaxis, patch.normal, indexes[i], fm.wSize, texs[i]))
      Normalize(texs[i]);

  inccs.resize(size);
  for (auto& incc : inccs)
		incc.resize(size);

  for (int i = 0; i < size; ++i) {
    inccs[i][i] = 0.0f;
    for (int j = i + 1; j < size; ++j) {
      if (!texs[i].empty() && !texs[j].empty())
        inccs[j][i] = inccs[i][j] = (!robust) ? 1.0f - Dot(texs[i], texs[j]) : RobustINCC(1.0f - Dot(texs[i], texs[j]));
      else
        inccs[j][i] = inccs[i][j] = 2.0f;
    }
  }
}

bool Optim::grabSafe(const int index, const int size, const Vec3f &center, const Vec3f &dx, const Vec3f &dy, const int level) const {
  const int margin = size / 2;

  const Vec3f tl = center - dx * margin - dy * margin;
  const Vec3f tr = center + dx * margin - dy * margin;

  const Vec3f bl = center - dx * margin + dy * margin;
  const Vec3f br = center + dx * margin + dy * margin;

  const float minx = min(tl[0], min(tr[0], min(bl[0], br[0])));
  const float maxx = max(tl[0], max(tr[0], max(bl[0], br[0])));
  const float miny = min(tl[1], min(tr[1], min(bl[1], br[1])));
  const float maxy = max(tl[1], max(tr[1], max(bl[1], br[1])));

  // 1 should be enough
  const int margin2 = 3;
  // ??? may need to change if we change interpolation method
  if (minx < margin2 || fm.ps.GetWidth(index,  level) - 1 - margin2 <= maxx ||
      miny < margin2 || fm.ps.GetHeight(index, level) - 1 - margin2 <= maxy)
    return false;
  return true;
}

// My own optimisaton
const float answers[] = {0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
float MyPow2(int x) {
  return answers[x + 4];
}

static float Log2 = log(2.0f);

bool Optim::grabTex(const Vec4f &coord, const Vec4f &pxaxis, const Vec4f &pyaxis, const Vec4f &pzaxis, 
										const int index, const int size, std::vector<float> &tex) const {
  tex.clear();

  Vec4f ray = fm.ps.photos[index].center - coord;
  unitize(ray);
  const float weight = max(0.0f, ray * pzaxis);

  //???????
  // if (weight < cos(fm.angleThreshold0))
  if (weight < cos(fm.angleThreshold1))
    return 1;

  const int margin = size / 2;

  Vec3f center = fm.ps.Project(index, coord,          fm.level);
  Vec3f dx     = fm.ps.Project(index, coord + pxaxis, fm.level) - center;
  Vec3f dy     = fm.ps.Project(index, coord + pyaxis, fm.level) - center;

  const float ratio = (norm(dx) + norm(dy)) / 2.0f;
  // int leveldif = (int)floor(log(ratio) / log(2.0f) + 0.5f);
  int leveldif = (int)floor(log(ratio) / Log2 + 0.5f);

  // Upper limit is 2
  leveldif = max(-fm.level, min(2, leveldif));

  // const float scale = pow(2.0f, (float)leveldif);

  const float scale    = MyPow2(leveldif);
  const int   newlevel = fm.level + leveldif;

  center /= scale;
  dx     /= scale;
  dy     /= scale;

  if (!grabSafe(index, size, center, dx, dy, newlevel))
    return true;

  Vec3f left = center - dx * margin - dy * margin;

  tex.resize(3 * size * size);
  float *texp = &tex[0] - 1;
  for (int y = 0; y < size; ++y) {
    Vec3f vftmp = left;
    left += dy;
    for (int x = 0; x < size; ++x) {
      Vec3f color = fm.ps.GetColor(index, vftmp[0], vftmp[1], newlevel);
      *(++texp) = color[0];
      *(++texp) = color[1];
      *(++texp) = color[2];
      vftmp += dx;
    }
  }

  return false;
}

double Optim::ComputeINCC(const Vec4f &coord, const Vec4f &normal, const std::vector<int> &indexes, const int id, const bool robust) {
  if ((int)indexes.size() < 2)
    return 2.0;

  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  GetPAxes(index, coord, normal, pxaxis, pyaxis);

  return computeINCC(coord, normal, indexes, pxaxis, pyaxis, id, robust);
}

double Optim::computeINCC(const Vec4f &coord, const Vec4f &normal, const std::vector<int> &indexes, 
													const Vec4f &pxaxis, const Vec4f &pyaxis, const int id, const bool robust) {
  if ((int)indexes.size() < 2)
    return 2.0;

  const int size = min(fm.tau, (int)indexes.size());
  vector<vector<float>> &texs = texsT[id];

  for (int i = 0; i < size; ++i) 
    if (!grabTex(coord, pxaxis, pyaxis, normal, indexes[i], fm.wSize, texs[i]))
      Normalize(texs[i]);

  if (texs[0].empty())
    return 2.0;

  double score = 0.0;

  // pure pairwise of reference based
#ifdef PMVS_PAIRNCC
  float totalweight = 0.0;
  for (int i = 0; i < size; ++i) {
    for (int j = i + 1; j < size; ++j) {
      if (!texs[i].empty() && !texs[j].empty()) {
        const float ftmp = weightsT[id][i] * weightsT[id][j];
        totalweight += ftmp;
        if (robust)
          score += robustincc(1.0 - dot(texs[i], texs[j])) * ftmp;
        else
          score += (1.0 - dot(texs[i], texs[j])) * ftmp;
      }
    }
  }

  if (totalweight == 0.0)
    score = 2.0;
  else
    score /= totalweight;
#else
  float totalweight = 0.0;
  for (int i = 1; i < size; ++i) {
    if (!texs[i].empty()) {
      totalweight += weightsT[id][i];
      if (robust)
        score += RobustINCC(1.0 - Dot(texs[0], texs[i])) * weightsT[id][i];
      else
        score += (1.0 - Dot(texs[0], texs[i])) * weightsT[id][i];
    }
  }
  if (totalweight == 0.0)
    score = 2.0;
  else
    score /= totalweight;
#endif

  return score;
}

void Optim::lfunc(double *p, double *hx, int m, int n, void *adata) {
  int iflag;
  one->func(n, m, p, hx, &iflag, adata);
}

void Optim::func(int m, int n, double *x, double *fvec, int *iflag, void *arg) {
  const int id = *((int *)arg);
  double xs[3] = {x[0], x[1], x[2]};

  for (int i = 0; i < m; ++i)
    fvec[i] = 2.0;

  Vec4f coord, normal;
  decode(coord, normal, xs, id);

  const int index = indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  GetPAxes(index, coord, normal, pxaxis, pyaxis);

  const int size = min(fm.tau, (int)indexesT[id].size());

  for (int i = 0; i < size; ++i) 
    if (!grabTex(coord, pxaxis, pyaxis, normal, indexesT[id][i], fm.wSize, texsT[id][i]))
      Normalize(texsT[id][i]);
  

  int count = -1;
  for (int i = 0; i < size; ++i) {
    for (int j = i + 1; j < size; ++j) {
      count++;
      if (texsT[id][i].empty() || texsT[id][j].empty())
        continue;

      fvec[count] = RobustINCC(1.0 - Dot(texsT[id][i], texsT[id][j]));
    }
  }
}

// Normalize only scale for each image
void Optim::Normalize(std::vector<std::vector<float>> &texs, const int size) {
  // compute average rgb
  Vec3f ave;
  int denom = 0;

  vector<Vec3f> rgbs;
  rgbs.resize(size);
  for (int i = 0; i < size; ++i) {
    if (texs[i].empty())
      continue;

    int count = 0;
    while (count < (int)texs[i].size()) {
      rgbs[i][0] += texs[i][count++];
      rgbs[i][1] += texs[i][count++];
      rgbs[i][2] += texs[i][count++];
    }
    rgbs[i] /= (int)texs[i].size() / 3;

    ave += rgbs[i];
    ++denom;
  }

  // overall average
  if (!denom)
    return;

  ave /= denom;

  // Scale all the colors
  for (int i = 0; i < size; ++i) {
    if (texs[i].empty())
      continue;
    int count = 0;
    // compute scale
    Vec3f scale;
    for (int j = 0; j < 3; ++j)
      if (rgbs[i][j] != 0.0f)
        scale[j] = ave[j] / rgbs[i][j];

    while (count < (int)texs[i].size()) {
      texs[i][count++] *= scale[0];
      texs[i][count++] *= scale[1];
      texs[i][count++] *= scale[2];
    }
  }
}

void Optim::Normalize(std::vector<float> &tex) {
  const int size  = (int)tex.size();
  const int size3 = size / 3;
  Vec3f ave;

  float *texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    ave[0] += *(++texp);
    ave[1] += *(++texp);
    ave[2] += *(++texp);
  }

  ave /= size3;

  float ave2 = 0.0;
  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    const float f0 = ave[0] - *(++texp);
    const float f1 = ave[1] - *(++texp);
    const float f2 = ave[2] - *(++texp);

    ave2 += f0 * f0 + f1 * f1 + f2 * f2;
  }

  ave2 = sqrt(ave2 / size);

  if (ave2 == 0.0f)
    ave2 = 1.0f;

  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    *(++texp) -= ave[0];
    *texp     /= ave2;
    *(++texp) -= ave[1];
    *texp     /= ave2;
    *(++texp) -= ave[2];
    *texp     /= ave2;
  }
}

float Optim::Dot(const std::vector<float> &tex0, const std::vector<float> &tex1) const {
#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like
  // begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i) 
    ans += tex0[i] * tex1[i];
  return ans / size;
#else
  const int size = (int)tex0.size();
  vector<float>::const_iterator i0 = tex0.begin();
  vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
    ans += (*i0) * (*i1) * template[i];
  }
  return ans;
#endif
}

float Optim::SSD(const std::vector<float> &tex0, const std::vector<float> &tex1) const {
  const float scale = 0.01;

#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like
  // begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i) 
    ans += fabs(tex0[i] - tex1[i]);
  return scale * ans / size;
#else
  const int size = (int)tex0.size();
  vector<float>::const_iterator i0 = tex0.begin();
  vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
    const float ftmp = fabs((*i0) - (*i1));
    // ans += (*i0) * (*i1) * template[i];
    ans += ftmp * template[i];
  }
  return scale * ans;
#endif
}

float Optim::GetUnit(const int index, const Vec4f &coord) const {
  const float fz   = norm(coord - fm.ps.photos[index].center);
  const float ftmp = ipScales[index];
  if (ftmp == 0.0)
    return 1.0;

  return 2.0 * fz * (0x0001 << fm.level) / ftmp;
}

// get x and y axis to collect textures given reference image and normal
void Optim::GetPAxes(const int index, const Vec4f &coord, const Vec4f &normal, Vec4f &pxaxis, Vec4f &pyaxis) const {
  // yasu changed here for fpmvs
  const float pscale = GetUnit(index, coord);

  Vec3f normal3(normal[0], normal[1], normal[2]);
  Vec3f yaxis3 = cross(normal3, xAxes[index]);
  unitize(yaxis3);
  Vec3f xaxis3 = cross(yaxis3, normal3);
  pxaxis[0] = xaxis3[0];
  pxaxis[1] = xaxis3[1];
  pxaxis[2] = xaxis3[2];
  pxaxis[3] = 0.0;
  pyaxis[0] = yaxis3[0];
  pyaxis[1] = yaxis3[1];
  pyaxis[2] = yaxis3[2];
  pyaxis[3] = 0.0;

  pxaxis *= pscale;
  pyaxis *= pscale;
  const float xdis = norm(fm.ps.Project(index, coord + pxaxis, fm.level) - fm.ps.Project(index, coord, fm.level));
  const float ydis = norm(fm.ps.Project(index, coord + pyaxis, fm.level) - fm.ps.Project(index, coord, fm.level));
  pxaxis /= xdis;
  pyaxis /= ydis;
}

void Optim::SetWeightsT(const ptch::Patch &patch, const int id) {
  ComputeUnits(patch, weightsT[id]);
  for (auto& w : weightsT[id])
    w = min(1.0f, weightsT[id][0] / w);
  weightsT[id][0] = 1.0f;
}
