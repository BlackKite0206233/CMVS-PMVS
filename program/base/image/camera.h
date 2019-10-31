#ifndef IMAGE_CAMERA_H
#define IMAGE_CAMERA_H

#include "../numeric/mat3.h"
#include "../numeric/mat4.h"
#include "../numeric/vec4.h"
#include <algorithm>
#include <climits>
#include <string>
#include <vector>

namespace img {

class Camera {
public:
  Camera(void);
  virtual ~Camera();

  // Update projection matrices from intrinsics and extrinsics
  void UpdateProjection(void);
  // Update all the camera related parameters
  void UpdateCamera(void);

  virtual void Init(const std::string cname, const int maxlevel);
  void Write(const std::string file);

  inline Vec3f Project(const Vec4f &coord, const int level) const;
  inline Vec3f Mult(const Vec4f &coord, const int level) const;

  static void SetProjection(const std::vector<float> &intrinsics, const std::vector<float> &extrinsics, std::vector<Vec4f> &projection, const int txtType);

  float GetScale(const Vec4f &coord, const int level) const;
  void GetPAxes(const Vec4f &coord, const Vec4f &normal, Vec4f &pxAxis, Vec4f &pyAxis, const int level = 0) const;

  void SetAxesScale(const float axesScale);

  static void Proj2Q(Mat4 &mat, double q[6]);
  static void Q2Proj(const double q[6], Mat4 &mat);
  static void SetProjectionSub(double params[], std::vector<Vec4f> &projection, const int level);

  float ComputeDistance(const Vec4f &point) const;
  float ComputeDepth(const Vec4f &point) const;
  float ComputeDepthDif(const Vec4f &rhs, const Vec4f &lhs) const;

  // Compute where the viewing ray passing through coord intersects
  // with the plane abcd.
  Vec4f Intersect(const Vec4f &coord, const Vec4f &abcd) const;
  void Intersect(const  Vec4f &coord, const Vec4f &abcd, Vec4f &cross, float &distance) const;
  // Computer a 3D coordinate that projects to a given image
  // coordinate. You can specify a different depth by the third
  // component of icoord.
  Vec4f Unproject(const Vec3f &iCoord, const int level) const;

  void SetK(Mat3f &K)   const;
  void SetRT(Mat4f &RT) const;

  void GetR(Mat3f &R) const;

  //----------------------------------------------------------------------
  // txt file name
  std::string cName;
  // Optical center
  Vec4f center;
  // Optical axis
  Vec4f oAxis;

  float ipScale;
  // 3x4 projection matrix
  std::vector<std::vector<Vec4f>> projection;
  Vec3f xAxis;
  Vec3f yAxis;
  Vec3f zAxis;

  // intrinsic and extrinsic camera parameters. Compact form.
  std::vector<float> intrinsics;
  std::vector<float> extrinsics;
  // camera parameter type
  int txtType;

protected:
  int maxLevel;

  float axesScale;

  Vec4f getOpticalCenter(void) const;
};

inline Vec3f Camera::Project(const Vec4f &coord, const int level) const {
  Vec3f vtmp;
  for (int i = 0; i < 3; ++i)
    vtmp[i] = projection[level][i] * coord;

  if (vtmp[2] <= 0.0) {
    vtmp[0] = -0xffff;
    vtmp[1] = -0xffff;
    vtmp[2] = -1.0f;
    return vtmp;
  } else
    vtmp /= vtmp[2];

  vtmp[0] = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), vtmp[0]));
  vtmp[1] = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), vtmp[1]));

  return vtmp;
};

inline Vec3f Camera::Mult(const Vec4f &coord, const int level) const {
  Vec3f vtmp;
  for (int i = 0; i < 3; ++i)
    vtmp[i] = projection[level][i] * coord;

  return vtmp;
};

template <class T> float ComputeEPD(const TMat3<T> &F, const TVec3<T> &p0, const TVec3<T> &p1) {
  TVec3<T> line = F * p1;
  const T ftmp = sqrt(line[0] * line[0] + line[1] * line[1]);
  if (ftmp == 0.0)
    return 0.0;

  line /= ftmp;
  return fabs(line * p0);
};

template <class T> void SetF(const img::Camera &lhs, const img::Camera &rhs, TMat3<T> &F, const int level = 0) {
  const TVec4<T> &p00 = lhs.projection[level][0];
  const TVec4<T> &p01 = lhs.projection[level][1];
  const TVec4<T> &p02 = lhs.projection[level][2];

  const TVec4<T> &p10 = rhs.projection[level][0];
  const TVec4<T> &p11 = rhs.projection[level][1];
  const TVec4<T> &p12 = rhs.projection[level][2];

  F[0][0] = det(TMat4<T>(p01, p02, p11, p12));
  F[0][1] = det(TMat4<T>(p01, p02, p12, p10));
  F[0][2] = det(TMat4<T>(p01, p02, p10, p11));

  F[1][0] = det(TMat4<T>(p02, p00, p11, p12));
  F[1][1] = det(TMat4<T>(p02, p00, p12, p10));
  F[1][2] = det(TMat4<T>(p02, p00, p10, p11));

  F[2][0] = det(TMat4<T>(p00, p01, p11, p12));
  F[2][1] = det(TMat4<T>(p00, p01, p12, p10));
  F[2][2] = det(TMat4<T>(p00, p01, p10, p11));
};

}; // namespace Image

#endif // CAMERA_H
