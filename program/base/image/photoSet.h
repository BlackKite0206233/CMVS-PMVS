#ifndef IMAGE_PHOTOSETS_H
#define IMAGE_PHOTOSETS_H

#include "photo.h"
#include <map>

namespace img {

class PhotoSet {
public:
  PhotoSet(void);
  virtual ~PhotoSet();

  void Init(const std::vector<int> &images, const std::string prefix, const int maxLevel, const int size, const bool alloc);

  // grabTex given 2D sampling information
  void GrabTex(const int index, const int level, const Vec2f &iCoord, const Vec2f &xAxis, const Vec2f &yAxis, 
			   std::vector<Vec3f> &tex, const int normalizef = 1) const;

  // grabTex given 3D sampling information
  void GrabTex(const int index, const int level, const Vec4f &coord,
               const Vec4f &pxAxis, const Vec4f &pyAxis, const Vec4f &pzAxis,
               std::vector<Vec3f> &tex, float &weight, const int normalizef = 1) const;

  void Write(const std::string outdir);
  void Free(void);
  void Free(const int level);

  void SetEdge(const float threshold);

  inline Vec3f Project(const int index, const Vec4f &coord, const int level) const;
  inline Vec3f Mult(const    int index, const Vec4f &coord, const int level) const;

  inline int GetWidth(const  int index, const int level) const;
  inline int GetHeight(const int index, const int level) const;

  inline Vec3f GetColor(const Vec4f &coord, const int index, const int level) const;
  inline Vec3f GetColor(const int index, const float fx, const float fy, const int level) const;
  inline Vec3f GetColor(const int index, const   int ix, const   int iy, const int level) const;

  inline bool GetMask(const Vec4f &coord, const int level) const;
  inline int  GetMask(const Vec4f &coord, const int index, const int level) const;
  inline int  GetMask(const int index, const float fx, const float fy, const int level) const;
  inline int  GetMask(const int index, const   int ix, const   int iy, const int level) const;

  inline int GetEdge(const Vec4f &coord, const int index, const int level) const;
  inline int GetEdge(const int index, const float fx, const float fy, const int level) const;
  inline int GetEdge(const int index, const   int ix, const   int iy, const int level) const;

  static float INCC(const std::vector<std::vector<Vec3f>> &texs, const std::vector<float> &weights);

  bool CheckAngles(const Vec4f &coord, const std::vector<int> &indexes, const float minAngle, const float maxAngle, const int tau) const;

  void GetMinMaxAngles(const Vec4f &coord, const std::vector<int> &indexes, float &minAngle, float &maxAngle) const;

  float ComputeDepth(const int index, const Vec4f &coord) const;

  // Take care of indexes
  std::vector<int>   images;
  std::vector<Photo> photos;

  int image2index(const int image) const;
  std::map<int, int> m_dict;

  // Number of cameras.
  int num;
  // Root directory
  std::string prefix;
  // maximum level
  int maxLevel;
  // Window size used to refine location
  int size;

  // getPAxes
  void GetPAxes(const int index, const Vec4f &coord, const Vec4f &normal, Vec4f &pxAxis, Vec4f &pyAxis) const;

  // pairwise distance based on optical center and viewing direction
  void SetDistances(void);
  std::vector<std::vector<float>> distances;

protected:
};

Vec3f PhotoSet::Project(const int index, const Vec4f &coord, const int level) const {
  return photos[index].Project(coord, level);
};

Vec3f PhotoSet::Mult(const int index, const Vec4f &coord, const int level) const {
  return photos[index].Mult(coord, level);
};

int PhotoSet::GetWidth(const int index, const int level) const {
  return photos[index].GetWidth(level);
};

int PhotoSet::GetHeight(const int index, const int level) const {
  return photos[index].GetHeight(level);
};

Vec3f PhotoSet::GetColor(const Vec4f &coord, const int index, const int level) const {
  return photos[index].GetColor(coord, level);
};

Vec3f PhotoSet::GetColor(const int index, const float fx, const float fy, const int level) const {
  return photos[index].Image::Image::GetColor(fx, fy, level);
};

Vec3f PhotoSet::GetColor(const int index, const int ix, const int iy, const int level) const {
  return photos[index].Image::Image::GetColor(ix, iy, level);
};

bool PhotoSet::GetMask(const Vec4f &coord, const int level) const {
  for (int index = 0; index < num; ++index)
    if (!GetMask(coord, index, level))
      return false;
  return true;
};

int PhotoSet::GetMask(const Vec4f &coord, const int index, const int level) const {
  return photos[index].GetMask(coord, level);
};

int PhotoSet::GetMask(const int index, const float fx, const float fy, const int level) const {
  return photos[index].Image::Image::GetMask(fx, fy, level);
};

int PhotoSet::GetMask(const int index, const int ix, const int iy, const int level) const {
  return photos[index].Image::Image::GetMask(ix, iy, level);
};

int PhotoSet::GetEdge(const Vec4f &coord, const int index, const int level) const {
  return photos[index].GetEdge(coord, level);
};

int PhotoSet::GetEdge(const int index, const float fx, const float fy, const int level) const {
  return photos[index].Image::Image::GetEdge(fx, fy, level);
};

int PhotoSet::GetEdge(const int index, const int ix, const int iy, const int level) const {
  return photos[index].Image::Image::GetEdge(ix, iy, level);
};

}; // namespace Image

#endif // IMAGE_PHOTOSETS_H
