#ifndef PMVS3_PATCH_H
#define PMVS3_PATCH_H

#include "../numeric/vec4.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>

namespace ptch {

class Patch {
public:
  Patch(void) {
    ncc     = -1.0;
    tImages = 0;
    fix     = 0;
    // dflag is initialized only once. if failed in one direction, we
    // never try that.
    dFlag   = 0;

    // All non-class member variables need to be initialized so that
    // they aren't just uninitialized memory.
    flag    = 0;
    id      = 0;
    dScale  = 0;
    aScale  = 0;
    tmp     = 0;
  }

  //----------------------------------------------------------------------
  // saved information
  // 3D coordinates of the center of the patch
  Vec4f coord;
  // patch outward normal vector
  Vec4f normal;

  // associated image ids. first image id is the reference one. images
  // can be non-targetting image.
  std::vector<int> images;
  std::vector<TVec2<int>> grids;

  // visible images. m_vimages must be targetting images.
  std::vector<int> vImages;
  std::vector<TVec2<int>> vGrids;

  //----------------------------------------------------------------------
  inline float Score(const float threshold) const {
    return std::max(0.0f, ncc - threshold) * (int)images.size();
  }
  inline float Score2(const float threshold) const {
    return std::max(0.0f, ncc - threshold) * tImages;
  }

  // average ncc
  float ncc;
  // number of targetting images in m_images
  int tImages;

  // flat for expansion
  // 0: not yet tested
  // 1: done
  int flag;

  // for directional flag
  unsigned char dFlag;

  // fixed patch or not
  char fix;

  // id number in m_ppatches
  int id;

  // scaling factor corresponding to one pixel difference
  float dScale;
  float aScale;

  float tmp;
};

typedef boost::shared_ptr<Patch> pPatch;

struct Spatchcmp {
  bool operator()(const pPatch& lhs, const pPatch& rhs) {
    if (lhs.get() < rhs.get())
      return true;
    else
      return false;
  }
};

std::istream &operator>>(std::istream &istr, ptch::Patch &rhs);
std::ostream &operator<<(std::ostream &ostr, const ptch::Patch &rhs);

}; // namespace Patch

#endif // PMVS3_PATCH_H
