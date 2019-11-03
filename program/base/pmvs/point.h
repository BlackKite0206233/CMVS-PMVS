#ifndef PMVS3_POINT_H
#define PMVS3_POINT_H

#include "../numeric/mat3.h"
#include "../numeric/vec4.h"

namespace PMVS3 {
class Point {
public:
  Point(void);
  virtual ~Point();

  Vec3f iCoord;
  float response;

  // 0: Harris
  // 1: DoG
  int type;

  // tempporary variable, used to store original imageid in initial match
  int iTmp;

  // 3D coordinate
  Vec4f coord;

  bool operator<(const Point &rhs) const {
    return response < rhs.response;
  }

  friend std::istream &operator>>(std::istream &istr, Point &rhs);
  friend std::ostream &operator<<(std::ostream &ostr, const Point &rhs);
};

bool SortPoint(const Point &a, const Point &b);

std::istream &operator>>(std::istream &istr, Point &rhs);
std::ostream &operator<<(std::ostream &ostr, const Point &rhs);
}; // namespace PMVS3

#endif // PMVS3_POINT_H
