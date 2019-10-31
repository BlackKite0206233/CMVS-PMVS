#ifndef PMVS3_SEED_H
#define PMVS3_SEED_H

#include "patch.h"
#include "point.h"
#include <boost/shared_ptr.hpp>
#include <vector>

#include "tinycthread.h"

namespace PMVS3 {
class FindMatch;
typedef boost::shared_ptr<Point> pPoint;

class Seed {
public:
  Seed(FindMatch &findMatch);
  virtual ~Seed(){};

  void Init(const std::vector<std::vector<Point>> &points);
  void Run(void);
  void Clear(void);

protected:
  void readPoints(const std::vector<std::vector<Point>> &points);
  bool canAdd(const int index, const int x, const int y);

  void initialMatch(const int index, const int id);
  void collectCells(const int index0, const int index1, const Point &p0, std::vector<Vec2i> &cells);

  void collectCandidates(const int index, const std::vector<int> &indexes, const Point &point, std::vector<pPoint> &vcp);

  bool initialMatchSub(const int index0, const int index1, const int id, ptch::Patch &patch);

  void unproject(const int index0, const int index1, const Point &p0, const Point &p1, Vec4f &coord) const;

  //----------------------------------------------------------------------
  FindMatch &fm;
  // points in a grid. For each index, grid
  std::vector<std::vector<std::vector<pPoint>>> pPoints;

  //----------------------------------------------------------------------
  // thread related
  //----------------------------------------------------------------------
  void initialMatchThread(void);
  static int initialMatchThreadTmp(void *arg);

  // Number of trials
  std::vector<int> sCounts;
  // Number of failures in the prep
  std::vector<int> fCounts0;
  // Number of failures in the post processing
  std::vector<int> fCounts1;
  // Number passes
  std::vector<int> pCounts;
};
}; // namespace PMVS3

#endif // PMVS3_SEED_H
