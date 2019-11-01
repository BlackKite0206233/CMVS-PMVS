#ifndef PMVS3_OPTIM_H
#define PMVS3_OPTIM_H

#include "patch.h"
#include <vector>

namespace PMVS3 {

class FindMatch;

class Optim {
public:
  Optim(FindMatch &findMatch);

  void Init(void);

  //-----------------------------------------------------------------
  // Image manipulation
  //-----------------------------------------------------------------
  void CollectImages(const int index, std::vector<int> &indexes) const;
  void AddImages(ptch::Patch &patch) const;
  void RemoveImagesEdge(ptch::Patch &patch) const;

  float GetUnit(const int index, const Vec4f &coord) const;

  void ComputeUnits(const ptch::Patch &patch, std::vector<int> &indexes, std::vector<float> &fineness, std::vector<Vec4f> &rays) const;
  void ComputeUnits(const ptch::Patch &patch, std::vector<float> &fineness) const;

  //-----------------------------------------------------------------
  // Optimization
  //-----------------------------------------------------------------

  bool PreProcess(ptch::Patch &patch, const int id, const int seed);
  void RefinePatch(ptch::Patch &patch, const int id, const int time);

  bool RefinePatchBFGS(ptch::Patch &patch, const int id, const int time, const int ncc);

  bool PostProcess(ptch::Patch &patch, const int id, const int seed);

  void SetRefImage(ptch::Patch &patch, const int id);

  bool Check(ptch::Patch &patch);

  std::vector<int> status;

protected:
  void filterImagesByAngle(ptch::Patch &patch);

  void sortImages(ptch::Patch &patch) const;
  void constraintImages(ptch::Patch       &patch, const float nccThreshold, const int id);
  void setRefConstraintImages(ptch::Patch &patch, const float nccThreshold, const int id);

  void setINCCs(const ptch::Patch &patch, std::vector<float>              &nccs, const std::vector<int> &indexes, const int id, const bool robust);

  void setINCCs(const ptch::Patch &patch, std::vector<std::vector<float>> &nccs, const std::vector<int> &indexes, const int id, const bool robust);

  bool grabTex(const Vec4f &coord, const Vec4f &pxaxis, const Vec4f &pyaxis, const Vec4f &pzaxis, 
							 const int index, const int size, std::vector<float> &tex) const;

  bool grabSafe(const int index, const int size, const Vec3f &center, const Vec3f &dx, const Vec3f &dy, const int level) const;

  /*
  double computeINCC(const Vec4f& coord, const Vec4f& normal,
                     const std::vector<int>& indexes, const int id,
                     const int robust);
  */
  double computeINCC(const Vec4f &coord, const Vec4f &normal, const std::vector<int> &indexes, 
										 const Vec4f &pxaxis, const Vec4f &pyaxis, const int id, const bool robust);

public:
  static void Normalize(std::vector<float> &tex);
  static void Normalize(std::vector<std::vector<float>> &texs, const int size);

  float Dot(const std::vector<float> &tex0, const std::vector<float> &tex1) const;
  float SSD(const std::vector<float> &tex0, const std::vector<float> &tex1) const;

protected:
  static void lfunc(double *p, double *hx, int m, int n, void *adata);
  void func(int m, int n, double *x, double *fvec, int *iflag, void *arg);

  // BFGS
  static double my_f(unsigned n, const double *x, double *grad, void *my_func_data);

  void encode(const Vec4f &coord, double *const vect, const int id) const;
  void encode(const Vec4f &coord, const Vec4f &normal, double *const vect, const int id) const;
  void decode(Vec4f &coord, Vec4f &normal, const double *const vect, const int id) const;
  void decode(Vec4f &coord, const double *const vect, const int id) const;

public:
  void SetWeightsT(const ptch::Patch &patch, const int id);

  double ComputeINCC(const Vec4f &coord, const Vec4f &normal, const std::vector<int> &indexes, const int id, const bool robust);
  void GetPAxes(const int index, const Vec4f &coord, const Vec4f &normal, Vec4f &pxaxis, Vec4f &pyaxis) const;

  static inline float RobustINCC(const float rhs) {
    return rhs / (1 + 3 * rhs);
  }

  static inline float UnrobustINCC(const float rhs) {
    return rhs / (1 - 3 * rhs);
  }

protected:
  void setAxesScales(void);

  static Optim *one;
  FindMatch &fm;

  //-----------------------------------------------------------------
  // Axes
  std::vector<Vec3f> xAxes;
  std::vector<Vec3f> yAxes;
  std::vector<Vec3f> zAxes;
  // Scales
  std::vector<float> ipScales;

  //-----------------------------------------------------------------
  // For threads
  std::vector<float> vect0T;
  std::vector<Vec4f> centersT;
  std::vector<Vec4f> raysT;
  std::vector<std::vector<int>> indexesT;
  std::vector<float> dScalesT;
  std::vector<float> aScalesT;

  // stores current parameters for derivative computation
  std::vector<Vec3f> paramsT;

  // Grabbed texture
  std::vector<std::vector<std::vector<float>>> texsT; // last is 7x7x3 patch
  // weights for refineDepthOrientationWeighed
  std::vector<std::vector<float>> weightsT;
  // Working array for levmar
  std::vector<std::vector<double>> worksT;
};
}; // namespace PMVS3

#endif // PMVS3_OPTIM_H
