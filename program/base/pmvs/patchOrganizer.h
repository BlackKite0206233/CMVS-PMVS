#ifndef PMVS3_PATCHORGANIZERS_H
#define PMVS3_PATCHORGANIZERS_H

#include "patch.h"
#include <queue>

namespace PMVS3 {

class FindMatch;

class P_compare {
public:
  bool operator()(const ptch::pPatch& lhs, const ptch::pPatch& rhs) const {
    return lhs->tmp < rhs->tmp;
  }
};

class PatchOrganizer {
public:
	PatchOrganizer(FindMatch &findMatch);

  void Init(void);
  void CollectPatches(const bool target = false);
  void CollectPatches(std::priority_queue<ptch::pPatch, std::vector<ptch::pPatch>, P_compare> &pqpatches);

  void CollectPatches(const int index, std::priority_queue<ptch::pPatch, std::vector<ptch::pPatch>, P_compare> &pqpatches);
  void CollectNonFixPatches(const int index, std::vector<ptch::pPatch> &ppatches);

  void WritePatches2(const std::string prefix, const bool bExportPLY, const bool bExportPatch, const bool bExportPSet);

  void WritePLY(const std::vector<ptch::pPatch> &patches, const std::string filename);
  void WritePLY(const std::vector<ptch::pPatch> &patches, const std::string filename, const std::vector<Vec3i> &colors);

  void ReadPatches(void);

  void ClearCounts(void);
  void ClearFlags(void);

  void SetGridsImages(ptch::Patch &patch, const std::vector<int> &images) const;
  void AddPatch(ptch::pPatch &ppatch);
  void RemovePatch(const ptch::pPatch &ppatch);
  void SetGrids(ptch::pPatch &ppatch) const;
  void SetGrids(ptch::Patch &patch) const;
  void SetVImagesVGrids(ptch::pPatch &ppatch);
  void SetVImagesVGrids(ptch::Patch &patch);
  void UpdateDepthMaps(ptch::pPatch &ppatch);

  bool IsVisible(const  ptch::Patch &patch, const int image, const int  ix, const int  iy, const float strict, const bool lock);
  bool IsVisible0(const ptch::Patch &patch, const int image,       int &ix,       int &iy, const float strict, const bool lock);

  void FindNeighbors(const ptch::Patch &patch, std::vector<ptch::pPatch> &neighbors,
                     const bool lock, const float scale = 1.0f, const int margin = 1, const bool skipvis = false);

  void SetScales(ptch::Patch &patch) const;

  float ComputeUnit(const ptch::Patch &patch) const;

  // change the contents of m_images from images to indexes
  void Image2Index(ptch::Patch &patch);
  // change the contents of m_images from indexes to images
  void Index2Image(ptch::Patch &patch);

  //----------------------------------------------------------------------
  // Widths of grids
  std::vector<int> gWidths;
  std::vector<int> gHeights;
  //----------------------------------------------------------------------
  // image, grid
  std::vector<std::vector<std::vector<ptch::pPatch>>> pGrids;
  // image, grid
  std::vector<std::vector<std::vector<ptch::pPatch>>> vpGrids;
  // Closest patch
  std::vector<std::vector<ptch::pPatch>> dpGrids;

  // all the patches in the current level of m_pgrids
  std::vector<ptch::pPatch> pPatches;

  // Check how many times patch optimization was performed for expansion
  std::vector<std::vector<unsigned char>> counts;

  static ptch::pPatch MAXDEPTH;
  static ptch::pPatch BACKGROUND;

protected:
  FindMatch &fm;
};
}; // namespace PMVS3

#endif // PMVS3_PATCHORGANIZERS_H
