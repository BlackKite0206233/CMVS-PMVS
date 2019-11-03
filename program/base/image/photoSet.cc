#include "photoSet.h"
#include <algorithm>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace img;

PhotoSet::PhotoSet(void) {}

PhotoSet::~PhotoSet() {}

void PhotoSet::Init(const std::vector<int> &images, const std::string prefix, const int maxlevel, const int size, const bool alloc) {
  this->images = images;
  num = (int)images.size();

  for (int i = 0; i < (int)images.size(); ++i)
    m_dict[images[i]] = i;

  this->prefix   = prefix;
  this->maxLevel = max(1, maxlevel);
  photos.resize(num);
  cerr << "Reading images: " << flush;
  for (int index = 0; index < num; ++index) {
    const int image = this->images[index];

    char test0[1024], test1[1024];
    char test2[1024], test3[1024];
    sprintf(test0, "%svisualize/%08d.ppm",  prefix.c_str(), image);
    sprintf(test1, "%svisualize/%08d.jpg",  prefix.c_str(), image);
    sprintf(test2, "%svisualize/%08d.png",  prefix.c_str(), image);
    sprintf(test3, "%svisualize/%08d.tiff", prefix.c_str(), image);
    if (ifstream(test0) || ifstream(test1)
#if defined(PMVS_HAVE_PNG)
        || ifstream(test2)
#endif
#if defined(PMVS_HAVE_TIFF)
        || ifstream(test3)
#endif
    ) {
      char name[1024], mname[1024], ename[1024], cname[1024];

      // Set name
      sprintf(name,  "%svisualize/%08d", prefix.c_str(), image);
      sprintf(mname, "%smasks/%08d",     prefix.c_str(), image);
      sprintf(ename, "%sedges/%08d",     prefix.c_str(), image);
      sprintf(cname, "%stxt/%08d.txt",   prefix.c_str(), image);

      photos[index].Init(name, mname, ename, cname, maxLevel);
      if (alloc)
        photos[index].Alloc();
      else
        photos[index].Alloc(1);
      cerr << '*' << flush;
    }
    // try 4 digits
    else {
      char name[1024], mname[1024], ename[1024], cname[1024];

      // Set name
      sprintf(name,  "%svisualize/%04d", prefix.c_str(), image);
      sprintf(mname, "%smasks/%04d",     prefix.c_str(), image);
      sprintf(ename, "%sedges/%04d",     prefix.c_str(), image);
      sprintf(cname, "%stxt/%04d.txt",   prefix.c_str(), image);

      photos[index].Init(name, mname, ename, cname, maxLevel);
      if (alloc)
        photos[index].Alloc();
      else
        photos[index].Alloc(1);
      cerr << '*' << flush;
    }

    /*
    const int image = m_images[index];
    char name[1024], mname[1024], ename[1024], cname[1024];

    // Set name
    sprintf(name, "%svisualize/%08d", prefix.c_str(), image);
    sprintf(mname, "%smasks/%08d", prefix.c_str(), image);
    sprintf(ename, "%sedges/%08d", prefix.c_str(), image);
    sprintf(cname, "%stxt/%08d.txt", prefix.c_str(), image);

    m_photos[index].init(name, mname, ename, cname, m_maxLevel);
    if (alloc)
      m_photos[index].alloc();
    else
      m_photos[index].alloc(1);
    cerr << '*' << flush;
    */
  }
  cerr << endl;
  const int margin = size / 2;
  this->size = 2 * margin + 1;
}

void PhotoSet::Free(void) {
  for (auto& p : photos)
    p.Free();
}

void PhotoSet::Free(const int level) {
  for (auto& p : photos)
		p.Free(level);
}

void PhotoSet::SetEdge(const float threshold) {
  for (auto& p : photos)
    p.SetEdge(threshold);
}

void PhotoSet::Write(const std::string outdir) {
  for (int index = 0; index < num; ++index) {
    const int image = images[index];
    char buffer[1024];
    sprintf(buffer, "%s%08d.txt", outdir.c_str(), image);

    photos[index].Write(buffer);
  }
}

// get x and y axis to collect textures given reference index and normal
void PhotoSet::GetPAxes(const int index, const Vec4f &coord, const Vec4f &normal, Vec4f &pxAxis, Vec4f &pyAxis) const {
  photos[index].GetPAxes(coord, normal, pxAxis, pyAxis);
}

void PhotoSet::GrabTex(const int index, const int level, const Vec2f &icoord, 
											 const Vec2f &xAxis, const Vec2f &yAxis, std::vector<Vec3f> &tex, const int normalizef) const {
  photos[index].GrabTex(level, icoord, xAxis, yAxis, size, tex, normalizef);
}

// grabTex given 3D sampling information
void PhotoSet::GrabTex(const int index, const int level, const Vec4f &coord,
                       const Vec4f &pxAxis, const Vec4f &pyAxis, const Vec4f &pzAxis, 
											 std::vector<Vec3f> &tex, float &weight, const int normalizef) const {
  photos[index].GrabTex(level, coord, pxAxis, pyAxis, pzAxis, size, tex, weight, normalizef);
}

float PhotoSet::INCC(const std::vector<std::vector<Vec3f>> &texs, const std::vector<float> &weights) {
  float incctmp = 0.0;
  float denom   = 0.0;
  for (int i = 0; i < (int)weights.size(); ++i) {
    if (texs[i].empty())
      continue;
    for (int j = i + 1; j < (int)weights.size(); ++j) {
      if (texs[j].empty())
        continue;

      const float weight = weights[i] * weights[j];
      const float ftmp   = Photo::Idot(texs[i], texs[j]);
      incctmp += ftmp * weight;
      denom   += weight;
    }
  }

  return (denom == 0.0) ? 2.0f : incctmp / denom;
}

void PhotoSet::GetMinMaxAngles(const Vec4f &coord, const std::vector<int> &indexes, float &minAngle, float &maxAngle) const {
  minAngle = M_PI;
  maxAngle = 0.0f;
  vector<Vec4f> rays;
  rays.resize((int)indexes.size());
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const int index = indexes[i];
    rays[i] = photos[index].center - coord;
    unitize(rays[i]);
  }

  for (int i = 0; i < (int)indexes.size(); ++i) {
    for (int j = i + 1; j < (int)indexes.size(); ++j) {
      const float dot   = max(-1.0f, min(1.0f, rays[i] * rays[j]));
      const float angle = acos(dot);
      minAngle = min(angle, minAngle);
      maxAngle = max(angle, maxAngle);
    }
  }
}

bool PhotoSet::CheckAngles(const Vec4f &coord, const std::vector<int> &indexes, const float minAngle, const float maxAngle, const int num) const {
  int count = 0;

  vector<Vec4f> rays;
  rays.resize((int)indexes.size());
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const int index = indexes[i];
    rays[i] = photos[index].center - coord;
    unitize(rays[i]);
  }

  for (int i = 0; i < (int)indexes.size(); ++i) {
    for (int j = i + 1; j < (int)indexes.size(); ++j) {
      const float dot   = max(-1.0f, min(1.0f, rays[i] * rays[j]));
      const float angle = acos(dot);
      if (minAngle < angle && angle < maxAngle)
        ++count;
    }
  }

  // if (count < num * (num - 1) / 2)
  return count < 1;
}

float PhotoSet::ComputeDepth(const int index, const Vec4f &coord) const {
  return photos[index].ComputeDepth(coord);
}

void PhotoSet::SetDistances(void) {
  distances.resize(num);
  float avedis = 0.0f;
  int denom    = 0;
  for (int i = 0; i < num; ++i) {
    distances[i].resize(num);
    for (int j = 0; j < num; ++j) {
      if (i == j)
        distances[i][j] = 0.0f;
      else {
        const float ftmp  = norm(photos[i].center - photos[j].center);
        distances[i][j] = ftmp;
        avedis += ftmp;
        denom++;
      }
    }
  }
  if (!denom)
    return;

  avedis /= denom;
  if (avedis == 0.0f) {
    cerr << "All the optical centers are identical..?" << endl;
    exit(1);
  }

  // plus angle difference
  for (int i = 0; i < num; ++i) {
    Vec4f ray0 = photos[i].oAxis;
    ray0[3] = 0.0f;
    for (int j = 0; j < num; ++j) {
      Vec4f ray1 = photos[j].oAxis;
      ray1[3] = 0.0f;

      distances[i][j] /= avedis;
      const float margin = cos(10.0f * M_PI / 180.0f);
      const float dis    = max(0.0f, 1.0f - ray0 * ray1 - margin);
      distances[i][j] += dis;
    }
  }
}

int PhotoSet::image2index(const int image) const {
  auto pos = m_dict.find(image);
  return (pos == m_dict.end()) ? -1 : pos->second;
}
