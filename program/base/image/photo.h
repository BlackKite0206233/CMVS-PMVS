#ifndef IMAGE_PHOTO_H
#define IMAGE_PHOTO_H

#include "../numeric/vec4.h"
#include "camera.h"
#include "image.h"

namespace img {

// Photo is an image with camera parameters
class Photo : public Image, public Camera {
public:
  Photo(void);
  virtual ~Photo();

  virtual void Init(const std::string name, const std::string mName, const std::string cName, const int maxLevel = 1);

  virtual void Init(const std::string name, const std::string mName, const std::string eName, const std::string cName, const int maxLevel = 1);

  // grabTex given 2D sampling information
  void GrabTex(const int level, const Vec2f &iCoord, const Vec2f &xAxis, const Vec2f &yAxis, 
			   const int size, std::vector<Vec3f> &tex, const bool normalizef = true) const;

  // grabTex given 3D sampling information
  void GrabTex(const int level, const Vec4f &coord, const Vec4f &pxAxis, const Vec4f &pyAxis, const Vec4f &pzAxis, 
               const int size, std::vector<Vec3f> &tex, float &weight, const bool normalizef = true) const;

  inline Vec3f GetColor(const float fx, const float fy, const int level) const;
  inline Vec3f GetColor(const Vec4f &coord, const int level) const;
  inline int GetMask(const Vec4f &coord, const int level) const;
  inline int GetEdge(const Vec4f &coord, const int level) const;

  static float Idot(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1);

  static void IdotC(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1, double *idc);

  static void Normalize(std::vector<Vec3f> &tex);

  static float SSD(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1);

protected:
};

Vec3f Photo::GetColor(const float fx, const float fy, const int level) const {
  return Image::GetColor(fx, fy, level);
};

Vec3f Photo::GetColor(const Vec4f &coord, const int level) const {
  const Vec3f icoord = Project(coord, level);
  return Image::GetColor(icoord[0], icoord[1], level);
};

int Photo::GetMask(const Vec4f &coord, const int level) const {
  if (masks[level].empty())
    return 1;

  const Vec3f icoord = Project(coord, level);
  return Image::GetMask(icoord[0], icoord[1], level);
};

int Photo::GetEdge(const Vec4f &coord, const int level) const {
  if (edges[level].empty())
    return 1;

  const Vec3f iCoord = Project(coord, level);

  if (iCoord[0] < 0 || widths[level] - 1 <= iCoord[0] || iCoord[1] < 0 || heights[level] - 1 <= iCoord[1])
    return 0;

  return Image::GetEdge(iCoord[0], iCoord[1], level);
};

}; // namespace Image

#endif // PHOTO_H
