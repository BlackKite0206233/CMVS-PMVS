#include "photo.h"
#include <fstream>

using namespace std;
using namespace img;

Photo::Photo(void) {}

Photo::~Photo() {}

void Photo::Init(const std::string name, const std::string mname, const std::string cname, const int maxLevel) {
  Image::Init(name, mname, maxLevel);
  Camera::Init(cname, maxLevel);
}

void Photo::Init(const std::string name, const std::string mname, const std::string ename, const std::string cname, const int maxLevel) {
  Image::Init(name, mname, ename, maxLevel);
  Camera::Init(cname, maxLevel);
}

float Photo::SSD(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1) {
  float ans = 0.0f;
  for (int i = 0; i < (int)tex0.size(); ++i)
    ans += norm2(tex0[i] - tex1[i]);

  // Make sure that the score is below 2.0f
  ans /= (int)tex0.size() * (255.0 * 255.0 * 3.0);

  return ans;
}

float Photo::Idot(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1) {
  if (tex0.empty() || tex1.empty()) {
    cerr << "Error in idot. Empty textures" << endl;
    exit(1);
  }
  float ans = 0.0;
  for (int i = 0; i < (int)tex0.size(); ++i) {
    ans += tex0[i] * tex1[i];
  }

  return 1.0f - ans / (3 * (int)tex0.size());
}

void Photo::IdotC(const std::vector<Vec3f> &tex0, const std::vector<Vec3f> &tex1, double *idc) {
  if (tex0.empty() || tex1.empty()) {
    cerr << "Error in idotC. Empty textures" << endl;
    exit(1);
  }
  idc[0] = 0.0;
  idc[1] = 0.0;
  idc[2] = 0.0;
  for (int i = 0; i < (int)tex0.size(); ++i) 
    for (int j = 0; j < 3; ++j)
      idc[j] += tex0[i][j] * tex1[i][j];

  for (int j = 0; j < 3; ++j)
    idc[j] = 1.0 - idc[j] / (int)tex0.size();
}

void Photo::Normalize(std::vector<Vec3f> &tex) {
  //----------------------------------------------------------------------
  // normalize average
  Vec3f ave;
  for (auto& t : tex)
    ave += t;
  ave /= (int)tex.size();

  for (auto& t : tex)
    t -= ave;
  //----------------------------------------------------------------------
  // compute variance
  float ave2 = 0.0f;
  for (auto& t : tex)
    ave2 += t * t;
  ave2 /= (int)tex.size() * 3;
  ave2  = sqrt(ave2);
  if (ave2 == 0.0f)
    ave2 = 1.0f;

  for (auto& t : tex)
    t /= ave2;
}

void Photo::GrabTex(const int level, const Vec2f &iCoord, const Vec2f &xAxis, const Vec2f &yAxis, 
										const int size, std::vector<Vec3f> &tex, const bool normalizef) const {
  const int margin = size / 2;

  // Check boundary condition
  const float maxx = iCoord[0] + size * fabs(xAxis[0]) + size * fabs(yAxis[0]);
  const float minx = iCoord[0] - size * fabs(xAxis[0]) - size * fabs(yAxis[0]);
  const float maxy = iCoord[1] + size * fabs(xAxis[1]) + size * fabs(yAxis[1]);
  const float miny = iCoord[1] - size * fabs(xAxis[1]) - size * fabs(yAxis[1]);

  tex.clear();
  if (minx < 0 || GetWidth(level)  - 1 <= maxx || 
	  miny < 0 || GetHeight(level) - 1 <= maxy)
    return;

  // tex.reserve(size * size);
  for (int y = -margin; y <= margin; ++y) {
    Vec2f v2ftmp = iCoord - margin * xAxis + y * yAxis;
    for (int x = -margin; x <= margin; ++x) {
      tex.push_back(Image::GetColor(v2ftmp[0], v2ftmp[1], level));
      v2ftmp += xAxis;
    }
  }

  if (normalizef)
    Normalize(tex);
}

void Photo::GrabTex(const int level, const Vec4f &coord,  const Vec4f &pxAxis, const Vec4f &pyAxis, const Vec4f &pzAxis, 
                    const int size, std::vector<Vec3f> &tex, float &weight, const bool normalizef) const {
  const int scale = 0x0001 << level;

  const Vec3f icoord3 = Project(coord, level);
  const Vec2f icoord(icoord3[0], icoord3[1]);

  const Vec3f xaxis3 = Project(coord + pxAxis * scale, level) - icoord3;
  const Vec2f xaxis(xaxis3[0], xaxis3[1]);

  const Vec3f yaxis3 = Project(coord + pyAxis * scale, level) - icoord3;
  const Vec2f yaxis(yaxis3[0], yaxis3[1]);

  GrabTex(level, icoord, xaxis, yaxis, size, tex, normalizef);

  Vec4f ray = center - coord;
  unitize(ray);
  weight = max(0.0f, pzAxis * ray);
}
