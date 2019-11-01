#ifndef PMVS3_DETECTOR_H
#define PMVS3_DETECTOR_H

#include "point.h"
#include <vector>
#include <set>

namespace PMVS3 {
class Detector {
public:
  static void SetGaussD(const float sigmaD, std::vector<float> &gaussD);
  static void SetGaussI(const float sigmaI, std::vector<float> &gaussI);
  virtual ~Detector() {}

protected:
  static float setThreshold(std::multiset<Point> &grid);
  bool isCloseBoundary(const int x, const int y, const int margin) const;
  int width;
  int height;
  std::vector<std::vector<Vec3f>> image;
  std::vector<std::vector<unsigned char>> mask;

public:
  template <class T>
  void ConvolveX(std::vector<std::vector<T>> &image, const std::vector<std::vector<unsigned char>> &mask, 
								 const std::vector<float> &filter, std::vector<std::vector<T>> &buffer) {
    const int width  = image[0].size();
    const int height = image.size();
    const int margin = ((int)filter.size()) / 2;

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        buffer[y][x] *= 0.0;

        if (!mask.empty() && !mask[y][x])
          continue;

        for (int j = 0; j < (int)filter.size(); ++j) {
          int xtmp = x + j - margin;
          if (xtmp < 0)
            xtmp = 0;
          else if (width <= xtmp)
            xtmp = width - 1;

          const int ytmp = y;

          if (!mask.empty() && !mask[ytmp][xtmp])
            continue;

          buffer[y][x] += filter[j] * image[ytmp][xtmp];
        }
      }
    }

    buffer.swap(image);
  }

  template <class T>
  void ConvolveY(std::vector<std::vector<T>> &image, const std::vector<std::vector<unsigned char>> &mask, 
				         const std::vector<float> &filter, std::vector<std::vector<T>> &buffer) {
    const int width  = image[0].size();
    const int height = image.size();
    const int margin = ((int)filter.size()) / 2;

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        buffer[y][x] *= 0.0;

        if (!mask.empty() && !mask[y][x])
          continue;

        for (int j = 0; j < (int)filter.size(); ++j) {
          const int xtmp = x;
          int ytmp = y + j - margin;
          if (ytmp < 0)
            ytmp = 0;
          else if (height <= ytmp)
            ytmp = height - 1;

          if (!mask.empty() && !mask[ytmp][xtmp])
            continue;

          buffer[y][x] += filter[j] * image[ytmp][xtmp];
        }
      }
    }

    buffer.swap(image);
  }

  template <class T>
  void ConvolveX(std::vector<std::vector<T>> &image, const std::vector<float> &filter, std::vector<std::vector<T>> &buffer) {
    const int width  = image[0].size();
    const int height = image.size();
    const int margin = ((int)filter.size()) / 2;

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        buffer[y][x] *= 0.0;

        for (int j = 0; j < (int)filter.size(); ++j) {
          const int xtmp = x + j - margin;
          const int ytmp = y;
          if (xtmp < 0 || width <= xtmp)
            continue;

          buffer[y][x] += filter[j] * image[ytmp][xtmp];
        }
      }
    }

    buffer.swap(image);
  }

  template <class T>
  void ConvolveY(std::vector<std::vector<T>> &image, const std::vector<float> &filter, std::vector<std::vector<T>> &buffer) {
    const int width  = image[0].size();
    const int height = image.size();
    const int margin = ((int)filter.size()) / 2;

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        buffer[y][x] *= 0.0;

        for (int j = 0; j < (int)filter.size(); ++j) {
          const int xtmp = x;
          const int ytmp = y + j - margin;
          if (ytmp < 0 || height <= ytmp)
            continue;

          buffer[y][x] += filter[j] * image[ytmp][xtmp];
        }
      }
    }

    buffer.swap(image);
  }
};
}; // namespace PMVS3

#endif // PMVS3_DETECTOR_H
