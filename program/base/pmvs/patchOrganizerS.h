#ifndef PMVS3_PATCHORGANIZERS_H
#define PMVS3_PATCHORGANIZERS_H

#include "patch.h"
#include <queue>

namespace PMVS3 {

class CfindMatch;

class P_compare {
public:
  bool operator()(const Patch::Ppatch &lhs, const Patch::Ppatch &rhs) const {
    return lhs->m_tmp < rhs->m_tmp;
  }
};

class CpatchOrganizerS {
public:
  CpatchOrganizerS(CfindMatch &findMatch);

  void init(void);
  void collectPatches(const bool target = false);
  void collectPatches(std::priority_queue<Patch::Ppatch, std::vector<Patch::Ppatch>, P_compare> &pqpatches);

  void collectPatches(const int index, std::priority_queue<Patch::Ppatch, std::vector<Patch::Ppatch>, P_compare> &pqpatches);
  void collectNonFixPatches(const int index, std::vector<Patch::Ppatch> &ppatches);

  void writePatches2(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet);

  void writePLY(const std::vector<Patch::Ppatch> &patches, const std::string filename);
  void writePLY(const std::vector<Patch::Ppatch> &patches, const std::string filename, const std::vector<Vec3i> &colors);

  void readPatches(void);

  void clearCounts(void);
  void clearFlags(void);

  void setGridsImages(Patch::Cpatch &patch, const std::vector<int> &images) const;
  void addPatch(Patch::Ppatch &ppatch);
  void removePatch(const Patch::Ppatch &ppatch);
  void setGrids(Patch::Ppatch &ppatch) const;
  void setGrids(Patch::Cpatch &patch) const;
  void setVImagesVGrids(Patch::Ppatch &ppatch);
  void setVImagesVGrids(Patch::Cpatch &patch);
  void updateDepthMaps(Patch::Ppatch &ppatch);

  bool isVisible(const  Patch::Cpatch &patch, const int image, const int &ix, const int &iy, const float strict, const bool lock);
  bool isVisible0(const Patch::Cpatch &patch, const int image,       int &ix,       int &iy, const float strict, const bool lock);

  void findNeighbors(const Patch::Cpatch &patch, std::vector<Patch::Ppatch> &neighbors, 
                     const bool lock, const float scale = 1.0f, const int margin = 1, const bool skipvis = false);

  void setScales(Patch::Cpatch &patch) const;

  float computeUnit(const Patch::Cpatch &patch) const;

  // change the contents of m_images from images to indexes
  void image2index(Patch::Cpatch &patch);
  // change the contents of m_images from indexes to images
  void index2image(Patch::Cpatch &patch);

  //----------------------------------------------------------------------
  // Widths of grids
  std::vector<int> m_gwidths;
  std::vector<int> m_gheights;
  //----------------------------------------------------------------------
  // image, grid
  std::vector<std::vector<std::vector<Patch::Ppatch>>> m_pgrids;
  // image, grid
  std::vector<std::vector<std::vector<Patch::Ppatch>>> m_vpgrids;
  // Closest patch
  std::vector<std::vector<Patch::Ppatch>> m_dpgrids;

  // all the patches in the current level of m_pgrids
  std::vector<Patch::Ppatch> m_ppatches;

  // Check how many times patch optimization was performed for expansion
  std::vector<std::vector<unsigned char>> m_counts;

  static Patch::Ppatch m_MAXDEPTH;
  static Patch::Ppatch m_BACKGROUND;

protected:
  CfindMatch &m_fm;
};
}; // namespace PMVS3

#endif // PMVS3_PATCHORGANIZERS_H
