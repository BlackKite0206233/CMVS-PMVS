#include "patch.h"
#include "../numeric/vec4.h"
#include <string>

using namespace std;

std::istream &ptch::operator>>(std::istream &istr, Patch &rhs) {
  string header;
  int itmp;
  istr >> header >> rhs.coord >> rhs.normal >> rhs.ncc >> rhs.dScale >> rhs.aScale;

  if (header == "PATCHA") {
    int type;
    Vec4f dir;
    istr >> type >> dir;
  }

  istr >> itmp;
  rhs.images.resize(itmp);
  for (auto& image : rhs.images)
    istr >> image;

  istr >> itmp;
  rhs.vImages.resize(itmp);
  for (auto& image : rhs.vImages)
    istr >> image;

  return istr;
}

std::ostream &ptch::operator<<(std::ostream &ostr, const Patch &rhs) {
  ostr << "PATCHS" << endl
       << rhs.coord << endl
       << rhs.normal << endl
       << rhs.ncc << ' ' << rhs.dScale << ' ' << rhs.aScale << endl
       << (int)rhs.images.size() << endl;
  for (const auto& image : rhs.images)
    ostr << image << ' ';
  ostr << endl;

  ostr << (int)rhs.vImages.size() << endl;
  for (const auto& image : rhs.vImages)
    ostr << image << ' ';
  ostr << endl;

  return ostr;
}
