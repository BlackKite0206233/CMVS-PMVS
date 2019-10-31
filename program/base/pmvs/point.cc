#include "point.h"
#include <iostream>

using namespace PMVS3;
using namespace std;

Point::Point(void) {
  response = -1.0;
  type = -1;
}

Point::~Point() {}

std::istream &PMVS3::operator>>(std::istream &istr, Point &rhs) {
  string header;
  char str[1024];
  istr >> str;
  header = string(str);
  istr >> rhs.iCoord[0] >> rhs.iCoord[1] >> rhs.response >> rhs.type;
  rhs.iCoord[2] = 1.0f;
  return istr;
}

std::ostream &PMVS3::operator<<(std::ostream &ostr, const Point &rhs) {
  ostr << "POINT0" << endl
       << rhs.iCoord[0] << ' ' << rhs.iCoord[1] << ' ' << rhs.response << ' ' << rhs.type;
  return ostr;
}

bool SortCpoint(const Point &a, const Point &b) {
  return a.response < b.response;
}
