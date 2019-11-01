#ifndef IMAGE_IMAGE_H
#define IMAGE_IMAGE_H

#include "../numeric/vec3.h"
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace img {

class Image {
public:
  Image(void);
  virtual ~Image();

  virtual void Init(const std::string name, const std::string mName, const int maxLevel = 1);
  virtual void Init(const std::string name, const std::string mName, const std::string ename, const int maxLevel = 1);

  void SIFT(const Vec3f &center, const Vec3f &xAxis, const Vec3f &yAxis, std::vector<float> &descriptor) const;

  void SIFT(const Vec3f &center, const Vec3f &xAxis, const Vec3f &yAxis, const int level, std::vector<float> &descriptor) const;

  void SetEdge(const float threshold);

  // access to image/masks
  inline Vec3f GetColor(const float fx, const float fy, const int level) const;
  inline Vec3f GetColor(const int   ix, const int   iy, const int level) const;

  inline void SetColor(const int ix, const int iy, const int level, const Vec3f &rgb);

  inline int GetMask(const float fx, const float fy, const int level) const;
  inline int GetMask(const int   ix, const int   iy, const int level) const;

  inline int GetEdge(const float fx, const float fy, const int level) const;
  inline int GetEdge(const int   ix, const int   iy, const int level) const;

  inline int GetWidth(const  int level = 0) const;
  inline int GetHeight(const int level = 0) const;
  // inline int getCWidth(const int beta, const int level = 0) const;
  // inline int getCHeight(const int beta, const int level = 0) const;

  inline const std::vector<unsigned char> &GetImage(const int level) const;
  inline const std::vector<unsigned char> &GetMask(const  int level) const;
  inline const std::vector<unsigned char> &GetEdge(const  int level) const;
  inline std::vector<unsigned char> &GetImage(const int level);
  inline std::vector<unsigned char> &GetMask(const  int level);
  inline std::vector<unsigned char> &GetEdge(const  int level);

  inline bool IsSafe(const Vec3f &iCoord, const int level) const;
  // Check if a mask image exists
  inline bool IsMask(void) const;
  // Check if an edge image exists
  inline bool IsEdge(void) const;

  // allocate and free memories, this function is also called when you call
  // getColor/getMask when the memory is not allocated
  void Alloc(const bool fast = false, const int filter = 0);
  // free memory
  void Free(void);
  // free memory below the specified level
  void Free(const int freeLevel);

  static bool ReadAnyImage(const std::string file,  std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static bool ReadPBMImage(const std::string file,  std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static bool WritePBMImage(const std::string file, std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static bool ReadPGMImage(const std::string file,  std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static bool WritePGMImage(const std::string file, const std::vector<unsigned char> &image, const int width, const int height);

  static bool ReadPPMImage(const std::string file, std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static bool WritePPMImage(const std::string file, const std::vector<unsigned char> &image, const int width, const int height);

  static bool ReadJpegImage(const std::string file, std::vector<unsigned char> &image, int &width, int &height, const bool fast);

  static void WriteJpegImage(const std::string file, const std::vector<unsigned char> &buffer, const int width, const int height, const bool flip = false);

  static float Hsdis(const float h0, const float s0, const float h1, const float s1);

  static void RGB2HSV(const float r, const float g, const float b, float &hr, float &sr, float &vr);
  static void RGB2HSV(const Vec3f &rgb, Vec3f &hsv);
  static void RGB2HSV(const Vec3f &rgb, float &hr, float &sr, float &vr);
  static void RGB2HSV(const float r, const float g, const float b, Vec3f &hsv);

  static void RGB2HS(const float r, const float g, const float b, float &h, float &s);
  static void RGB2HS(const Vec3f &rgb, float &h, float &s);
  static void RGB2HS(const Vec3f &rgb, Vec2f &hs);
  static void RGB2HS(const float r, const float g, const float b, Vec2f &hs);

  static void Gray2RGB(const float gray, float &r, float &g, float &b);

  // Some low-level image processing
  // 2D convolution with twice 1D gaussian convolution.

  // Create a 1d gaussian filter based on sigma
  static void CreateFilter(const float sigma, std::vector<float> &filter);

  static void FilterG(const std::vector<float> &filter, std::vector<std::vector<float>> &data);
  static void FilterG(const std::vector<float> &filter, std::vector<std::vector<float>> &data, std::vector<std::vector<float>> &buffer);

  static void FilterG(const std::vector<float> &filter, const int width, const int height, std::vector<float> &data);

  static void FilterG(const std::vector<float> &filter, const int width, const int height, std::vector<float> &data, std::vector<float> &buffer);

  // non maximum surpression
  static void NMS(std::vector<std::vector<float>> &data);
  static void NMS(std::vector<std::vector<float>> &data, std::vector<std::vector<float>> &buffer);

  /*
  // ThreshMode
  // Used to filter out outliers. Very general algorithm.
  static void setInOut(const std::vector<std::vector<float> >& data,
  std::vector<int>& inout, const float sigma = 1.0f, const int specular = 0);

  // Used to filter out outliers. RGB version for specular highlights.
  static void setInOut(const std::vector<Vec3f>& rgbs, std::vector<int>& inout,
                       const float sigma = 1.0f, const int specular = 0);
  

  // Used to filter out outliers. HSV version for specular highlights.
  static void setInOutHSV(const std::vector<Vec3f>& hsvs, std::vector<int>&
  inout, const float sigma = 1.0f, const int specular = 0);
  */
protected:
  //----------------------------------------------------------------------
  // member functions
  //----------------------------------------------------------------------
  // complete the name of an image file
  static void completeName(const std::string &lhs, std::string &rhs, const bool color);

  // build image pyramids
  void buildImageMaskEdge(const int filter);
  void buildImage(const int filter);
  void buildMask(void);
  void buildEdge(void);

#ifdef FURUKAWA_IMAGE_GAMMA
  void decodeGamma(void);
#endif

  //----------------------------------------------------------------------
  // Variables updated at every alloc/free
  //----------------------------------------------------------------------
  // 0: nothing allocated
  // 1: width/height allocated
  // 2: memory allocated
  int alloc;
  // a pyramid of images
  std::vector<std::vector<unsigned char>> images;
  // a pyramid of masks
  std::vector<std::vector<unsigned char>> masks;
  // a pyramid of images specifying regions with edges(texture)
  std::vector<std::vector<unsigned char>> edges;

  // width of an image in each level
  std::vector<int> widths;
  // height of an image in each level
  std::vector<int> heights;

#ifdef FURUKAWA_IMAGE_GAMMA
  // For gamma decoded images
  std::vector<std::vector<float>> m_dimages;
#endif

  //----------------------------------------------------------------------
  // Variables keep fixed
  //----------------------------------------------------------------------
  // a name of an image
  std::string name;
  // a name of a mask image
  std::string mName;
  // a name of an image specifying regions with edges(texture)
  std::string eName;

  // number of levels
  int maxLevel;
};

inline bool Image::IsSafe(const Vec3f &iCoord, const int level) const {
#ifdef FURUKAWA_IMAGE_BICUBIC
  if (icoord[0] < 1.0 || widths[level] - 3 < icoord[0] || icoord[1] < 1.0 ||
      heights[level] - 3 < icoord[1])
    return 0;
  else
    return 1;
#else
  return !(iCoord[0] < 0.0 || widths[level] - 2 < iCoord[0] || iCoord[1] < 0.0 || heights[level] - 2 < iCoord[1]);
#endif
};

// Check if a mask image exists
inline bool Image::IsMask(void) const {
  return !masks[0].empty();
};

// Check if an edge image exists
inline bool Image::IsEdge(void) const {
  return !edges[0].empty();
};

inline const std::vector<unsigned char> & Image::GetImage(const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

#ifdef FURUKAWA_IMAGE_GAMMA
  std::cerr << "Cannot do getImage in the gamma correction mode." << std::endl;
  exit(1);
#endif

  return images[level];
};

inline const std::vector<unsigned char> & Image::GetMask(const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
  return masks[level];
};

inline const std::vector<unsigned char> & Image::GetEdge(const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
  return edges[level];
};

inline std::vector<unsigned char> &Image::GetImage(const int level) {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

#ifdef FURUKAWA_IMAGE_GAMMA
  std::cerr << "Cannot do getImage in the gamma correction mode." << std::endl;
  exit(1);
#endif

  return images[level];
};

inline std::vector<unsigned char> &Image::GetMask(const int level) {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
  return masks[level];
};

inline std::vector<unsigned char> &Image::GetEdge(const int level) {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
  return edges[level];
};

int Image::GetWidth(const int level) const {
  if (alloc == 0) {
    std::cerr << "First allocate (getWidth)" << std::endl;
    exit(1);
  }
  return widths[level];
};

int Image::GetHeight(const int level) const {
  if (alloc == 0) {
    // alloc(1);
    std::cerr << "First allocate (getHeight)" << std::endl;
    exit(1);
  }
  return heights[level];
};

/*
int Image::getCWidth(const int beta, const int level) const{
  if (alloc == 0) {
    //alloc(1);
    std::cerr << "First allocate (getCWidth)" << std::endl;
    exit (1);
  }
  return (widths[level] - 1) / beta + 1;
};
*/
/*
int Image::getCHeight(const int beta, const int level) const{
  if (alloc == 0) {
    //alloc(1);
    std::cerr << "First allocate (getCHeight)" << std::endl;
    exit (1);
  }
  return (heights[level] - 1) / beta + 1;
};
*/

Vec3f Image::GetColor(const float x, const float y, const int level) const {
#ifdef FURUKAWA_DEBUG
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
#endif

#ifdef FURUKAWA_IMAGE_BICUBIC
  const int x1 = (int)floor(x);
  const int y1 = (int)floor(y);
  const float p = x - x1;
  const float q = y - y1;

  float f = 1 + p;
  const float wx0 = (((-1) * f + 5) * f - 8) * f + 4;
  f = 2 - p;
  const float wx3 = (((-1) * f + 5) * f - 8) * f + 4;
  f = p;
  const float wx1 = (((1) * f - 2) * f) * f + 1;
  f = 1 - p;
  const float wx2 = (((1) * f - 2) * f) * f + 1;

  f = 1 + q;
  const float wy0 = (((-1) * f + 5) * f - 8) * f + 4;
  f = 2 - q;
  const float wy3 = (((-1) * f + 5) * f - 8) * f + 4;
  f = q;
  const float wy1 = (((1) * f - 2) * f) * f + 1;
  f = 1 - q;
  const float wy2 = (((1) * f - 2) * f) * f + 1;

  const int offset = widths[level] * 3;
  const int index0 = ((y1 - 1) * widths[level] + x1 - 1) * 3;
  const int index1 = index0 + offset;
  const int index2 = index1 + offset;
  const int index3 = index2 + offset;

#ifdef FURUKAWA_IMAGE_GAMMA
  const float &r00 = m_dimages[level][index0];
  const float &g00 = m_dimages[level][index0 + 1];
  const float &b00 = m_dimages[level][index0 + 2];
  const float &r01 = m_dimages[level][index0 + 3];
  const float &g01 = m_dimages[level][index0 + 4];
  const float &b01 = m_dimages[level][index0 + 5];
  const float &r02 = m_dimages[level][index0 + 6];
  const float &g02 = m_dimages[level][index0 + 7];
  const float &b02 = m_dimages[level][index0 + 8];
  const float &r03 = m_dimages[level][index0 + 9];
  const float &g03 = m_dimages[level][index0 + 10];
  const float &b03 = m_dimages[level][index0 + 11];

  const float &r10 = m_dimages[level][index1];
  const float &g10 = m_dimages[level][index1 + 1];
  const float &b10 = m_dimages[level][index1 + 2];
  const float &r11 = m_dimages[level][index1 + 3];
  const float &g11 = m_dimages[level][index1 + 4];
  const float &b11 = m_dimages[level][index1 + 5];
  const float &r12 = m_dimages[level][index1 + 6];
  const float &g12 = m_dimages[level][index1 + 7];
  const float &b12 = m_dimages[level][index1 + 8];
  const float &r13 = m_dimages[level][index1 + 9];
  const float &g13 = m_dimages[level][index1 + 10];
  const float &b13 = m_dimages[level][index1 + 11];

  const float &r20 = m_dimages[level][index2];
  const float &g20 = m_dimages[level][index2 + 1];
  const float &b20 = m_dimages[level][index2 + 2];
  const float &r21 = m_dimages[level][index2 + 3];
  const float &g21 = m_dimages[level][index2 + 4];
  const float &b21 = m_dimages[level][index2 + 5];
  const float &r22 = m_dimages[level][index2 + 6];
  const float &g22 = m_dimages[level][index2 + 7];
  const float &b22 = m_dimages[level][index2 + 8];
  const float &r23 = m_dimages[level][index2 + 9];
  const float &g23 = m_dimages[level][index2 + 10];
  const float &b23 = m_dimages[level][index2 + 11];

  const float &r30 = m_dimages[level][index3];
  const float &g30 = m_dimages[level][index3 + 1];
  const float &b30 = m_dimages[level][index3 + 2];
  const float &r31 = m_dimages[level][index3 + 3];
  const float &g31 = m_dimages[level][index3 + 4];
  const float &b31 = m_dimages[level][index3 + 5];
  const float &r32 = m_dimages[level][index3 + 6];
  const float &g32 = m_dimages[level][index3 + 7];
  const float &b32 = m_dimages[level][index3 + 8];
  const float &r33 = m_dimages[level][index3 + 9];
  const float &g33 = m_dimages[level][index3 + 10];
  const float &b33 = m_dimages[level][index3 + 11];
#else
  const unsigned char &r00 = images[level][index0];
  const unsigned char &g00 = images[level][index0 + 1];
  const unsigned char &b00 = images[level][index0 + 2];
  const unsigned char &r01 = images[level][index0 + 3];
  const unsigned char &g01 = images[level][index0 + 4];
  const unsigned char &b01 = images[level][index0 + 5];
  const unsigned char &r02 = images[level][index0 + 6];
  const unsigned char &g02 = images[level][index0 + 7];
  const unsigned char &b02 = images[level][index0 + 8];
  const unsigned char &r03 = images[level][index0 + 9];
  const unsigned char &g03 = images[level][index0 + 10];
  const unsigned char &b03 = images[level][index0 + 11];

  const unsigned char &r10 = images[level][index1];
  const unsigned char &g10 = images[level][index1 + 1];
  const unsigned char &b10 = images[level][index1 + 2];
  const unsigned char &r11 = images[level][index1 + 3];
  const unsigned char &g11 = images[level][index1 + 4];
  const unsigned char &b11 = images[level][index1 + 5];
  const unsigned char &r12 = images[level][index1 + 6];
  const unsigned char &g12 = images[level][index1 + 7];
  const unsigned char &b12 = images[level][index1 + 8];
  const unsigned char &r13 = images[level][index1 + 9];
  const unsigned char &g13 = images[level][index1 + 10];
  const unsigned char &b13 = images[level][index1 + 11];

  const unsigned char &r20 = images[level][index2];
  const unsigned char &g20 = images[level][index2 + 1];
  const unsigned char &b20 = images[level][index2 + 2];
  const unsigned char &r21 = images[level][index2 + 3];
  const unsigned char &g21 = images[level][index2 + 4];
  const unsigned char &b21 = images[level][index2 + 5];
  const unsigned char &r22 = images[level][index2 + 6];
  const unsigned char &g22 = images[level][index2 + 7];
  const unsigned char &b22 = images[level][index2 + 8];
  const unsigned char &r23 = images[level][index2 + 9];
  const unsigned char &g23 = images[level][index2 + 10];
  const unsigned char &b23 = images[level][index2 + 11];

  const unsigned char &r30 = images[level][index3];
  const unsigned char &g30 = images[level][index3 + 1];
  const unsigned char &b30 = images[level][index3 + 2];
  const unsigned char &r31 = images[level][index3 + 3];
  const unsigned char &g31 = images[level][index3 + 4];
  const unsigned char &b31 = images[level][index3 + 5];
  const unsigned char &r32 = images[level][index3 + 6];
  const unsigned char &g32 = images[level][index3 + 7];
  const unsigned char &b32 = images[level][index3 + 8];
  const unsigned char &r33 = images[level][index3 + 9];
  const unsigned char &g33 = images[level][index3 + 10];
  const unsigned char &b33 = images[level][index3 + 11];
#endif
  // separate x and y
  const float row0[3] = {wx0 * r00 + wx1 * r01 + wx2 * r02 + wx3 * r03,
                         wx0 * g00 + wx1 * g01 + wx2 * g02 + wx3 * g03,
                         wx0 * b00 + wx1 * b01 + wx2 * b02 + wx3 * b03};
  const float row1[3] = {wx0 * r10 + wx1 * r11 + wx2 * r12 + wx3 * r13,
                         wx0 * g10 + wx1 * g11 + wx2 * g12 + wx3 * g13,
                         wx0 * b10 + wx1 * b11 + wx2 * b12 + wx3 * b13};
  const float row2[3] = {wx0 * r20 + wx1 * r21 + wx2 * r22 + wx3 * r23,
                         wx0 * g20 + wx1 * g21 + wx2 * g22 + wx3 * g23,
                         wx0 * b20 + wx1 * b21 + wx2 * b22 + wx3 * b23};
  const float row3[3] = {wx0 * r30 + wx1 * r31 + wx2 * r32 + wx3 * r33,
                         wx0 * g30 + wx1 * g31 + wx2 * g32 + wx3 * g33,
                         wx0 * b30 + wx1 * b31 + wx2 * b32 + wx3 * b33};

  float r = wy0 * row0[0] + wy1 * row1[0] + wy2 * row2[0] + wy3 * row3[0];
  float g = wy0 * row0[1] + wy1 * row1[1] + wy2 * row2[1] + wy3 * row3[1];
  float b = wy0 * row0[2] + wy1 * row1[2] + wy2 * row2[2] + wy3 * row3[2];

  return Vec3f(r, g, b);
#else
  // Bilinear case
  const int lx    = static_cast<int>(x);
  const int ly    = static_cast<int>(y);
  const int index = 3 * (ly * widths[level] + lx);

  const float dx1 =    x - lx;
  const float dx0 = 1.0f - dx1;
  const float dy1 =    y - ly;
  const float dy0 = 1.0f - dy1;

  const float f00  = dx0 * dy0;
  const float f01  = dx0 * dy1;
  const float f10  = dx1 * dy0;
  const float f11  = dx1 * dy1;
  const int index2 = index + 3 * widths[level];

#ifdef FURUKAWA_IMAGE_GAMMA
  const float *fp0 = &m_dimages[level][index] - 1;
  const float *fp1 = &m_dimages[level][index2] - 1;
  float r = 0.0f;
  float g = 0.0f;
  float b = 0.0f;
  r += *(++fp0) * f00 + *(++fp1) * f01;
  g += *(++fp0) * f00 + *(++fp1) * f01;
  b += *(++fp0) * f00 + *(++fp1) * f01;
  r += *(++fp0) * f10 + *(++fp1) * f11;
  g += *(++fp0) * f10 + *(++fp1) * f11;
  b += *(++fp0) * f10 + *(++fp1) * f11;
  return Vec3f(r, g, b);
  /*
  return Vec3f(m_dimages[level][index] * f00 + m_dimages[level][index + 3] * f10
  + m_dimages[level][index2] * f01 + m_dimages[level][index2 + 3] * f11,
               m_dimages[level][index + 1] * f00 + m_dimages[level][index + 4] *
  f10 + m_dimages[level][index2 + 1] * f01 + m_dimages[level][index2 + 4] * f11,
               m_dimages[level][index + 2] * f00 + m_dimages[level][index + 5] *
  f10 + m_dimages[level][index2 + 2] * f01 + m_dimages[level][index2 + 5] *
  f11);
  */
#else
  const unsigned char *ucp0 = &images[level][index]  - 1;
  const unsigned char *ucp1 = &images[level][index2] - 1;
  float r = 0.0f;
  float g = 0.0f;
  float b = 0.0f;
  r += *(++ucp0) * f00 + *(++ucp1) * f01;
  g += *(++ucp0) * f00 + *(++ucp1) * f01;
  b += *(++ucp0) * f00 + *(++ucp1) * f01;
  r += *(++ucp0) * f10 + *(++ucp1) * f11;
  g += *(++ucp0) * f10 + *(++ucp1) * f11;
  b += *(++ucp0) * f10 + *(++ucp1) * f11;
  return Vec3f(r, g, b);
  /*
  return Vec3f(images[level][index] * f00 + images[level][index + 3] * f10 +
               images[level][index2] * f01 + images[level][index2 + 3] *
  f11,
               

               images[level][index + 1] * f00 + images[level][index + 4] *
  f10 + images[level][index2 + 1] * f01 + images[level][index2 + 4] * f11,
               

               images[level][index + 2] * f00 + images[level][index + 5] *
  f10 + images[level][index2 + 2] * f01 + images[level][index2 + 5] * f11);
  */
#endif
  /*
  const int lx = (int)floor(x);    const int ux = lx + 1;
  const int ly = (int)floor(y);    const int uy = ly + 1;

  const Vec3f vc[2][2] = {{getColor(lx, ly, level), getColor(lx, uy, level)},
                          {getColor(ux, ly, level), getColor(ux, uy, level)}};
  

  const Vec3f color = vc[0][0] * (ux - x) * (uy - y) + vc[0][1] * (ux - x) * (y
  - ly) + vc[1][0] * (x - lx) * (uy - y) + vc[1][1] * (x - lx) * (y - ly);
  return color;
  */
#endif
};

void Image::SetColor(const int ix, const int iy, const int level, const Vec3f &rgb) {
#ifdef FURUKAWA_DEBUG
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
#endif
  const int index = (iy * widths[level] + ix) * 3;

#ifdef FURUKAWA_IMAGE_GAMMA
  m_dimages[level][index] = rgb[0];
  m_dimages[level][index + 1] = rgb[1];
  m_dimages[level][index + 2] = rgb[2];
#else
  images[level][index    ] = (unsigned char)floor(rgb[0] + 0.5f);
  images[level][index + 1] = (unsigned char)floor(rgb[1] + 0.5f);
  images[level][index + 2] = (unsigned char)floor(rgb[2] + 0.5f);
#endif
};

Vec3f Image::GetColor(const int ix, const int iy, const int level) const {
#ifdef FURUKAWA_DEBUG
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }
#endif
  const int index = (iy * widths[level] + ix) * 3;

#ifdef FURUKAWA_IMAGE_GAMMA
  return Vec3f(m_dimages[level][index], m_dimages[level][index + 1],
               m_dimages[level][index + 2]);
#else
  return Vec3f(images[level][index], images[level][index + 1], images[level][index + 2]);
#endif
};

int Image::GetMask(const float fx, const float fy, const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

  if (masks[level].empty())
    return 1;

  const int ix = (int)floor(fx + 0.5f);
  const int iy = (int)floor(fy + 0.5f);
  return GetMask(ix, iy, level);
};

int Image::GetMask(const int ix, const int iy, const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

  if (masks[level].empty())
    return 1;

  if (ix < 0 || widths[level] <= ix || iy < 0 || heights[level] <= iy)
    return 1;

  const int index = iy * widths[level] + ix;
  return masks[level][index];
};

int Image::GetEdge(const float fx, const float fy, const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

  if (edges[level].empty())
    return 1;
  const int ix = (int)floor(fx + 0.5f);
  const int iy = (int)floor(fy + 0.5f);
  return GetEdge(ix, iy, level);
};

int Image::GetEdge(const int ix, const int iy, const int level) const {
  if (alloc != 2) {
    std::cerr << "First allocate" << std::endl;
    exit(1);
  }

  if (edges[level].empty())
    return 1;

  if (ix < 0 || widths[level] <= ix || iy < 0 || heights[level] <= iy)
    return 1;

  const int index = iy * widths[level] + ix;
  return edges[level][index];
};

}; // namespace Image

#endif // IMAGE_H
