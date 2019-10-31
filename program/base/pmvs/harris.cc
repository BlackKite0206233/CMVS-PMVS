#include "harris.h"
#include <algorithm>

using namespace PMVS3;
using namespace std;

void Harris::init(const std::vector<unsigned char> &image, const std::vector<unsigned char> &mask, const std::vector<unsigned char> &edge) {
  this->image.clear();
  this->image.resize(height);
  int count = 0;
  for (int y = 0; y < height; ++y) {
    this->image[y].resize(width);
    for (int x = 0; x < width; ++x) {
	  this->image[y][x][0] = ((int)image[count++]) / 255.0f;
	  this->image[y][x][1] = ((int)image[count++]) / 255.0f;
	  this->image[y][x][2] = ((int)image[count++]) / 255.0f;
    }
  }

  this->mask.clear();
  if (!mask.empty() || !edge.empty()) {
    this->mask.resize(height);
    count = 0;
    for (int y = 0; y < height; ++y) {
      this->mask[y].resize(width);
      for (int x = 0; x < width; ++x) {
        if (mask.empty())
          this->mask[y][x] = edge[count++];
        else if (edge.empty())
          this->mask[y][x] = mask[count++];
        else {
          if (mask[count] && edge[count])
            this->mask[y][x] = (unsigned char)255;
          else
            this->mask[y][x] = 0;
          count++;
        }
      }
    }
  }

  SetGaussD(sigmaD, gaussD);
  SetGaussI(sigmaI, gaussI);
}

void Harris::setDerivatives(void) {
  // set dIdx, dIdy
  preprocess();

  // now set dIdxdIdx, dIdydIdy, dIdxDIdy
  preprocess2();
}

void Harris::preprocess2(void) {
  dIdxdIdx.clear();
  dIdydIdy.clear();
  dIdxdIdy.clear();
  dIdxdIdx.resize(height);
  dIdydIdy.resize(height);
  dIdxdIdy.resize(height);
  for (int y = 0; y < height; ++y) {
    dIdxdIdx[y].resize(width);
    dIdydIdy[y].resize(width);
    dIdxdIdy[y].resize(width);
    for (int x = 0; x < width; ++x) {
      dIdxdIdx[y][x] = dIdydIdy[y][x] = dIdxdIdy[y][x] = 0.0;
      if (!mask.empty() && !mask[y][x])
        continue;

      dIdxdIdx[y][x] += dIdx[y][x] * dIdx[y][x];
      dIdydIdy[y][x] += dIdy[y][x] * dIdy[y][x];
      dIdxdIdy[y][x] += dIdx[y][x] * dIdy[y][x];
    }
  }

  {
    vector<vector<Vec3f>>().swap(dIdx);
    vector<vector<Vec3f>>().swap(dIdy);
  }

  //----------------------------------------------------------------------
  // blur
  vector<vector<float>> vvftmp;
  vvftmp.resize(height);
  for (int y = 0; y < height; ++y) {
    vvftmp[y].resize(width);
    for (int x = 0; x < width; ++x)
      vvftmp[y][x] = 0.0;
  }

  //----------------------------------------------------------------------
  // dIdxdIdx
  ConvolveX(dIdxdIdx, mask, gaussI, vvftmp);
  ConvolveY(dIdxdIdx, mask, gaussI, vvftmp);

  //----------------------------------------------------------------------
  // dIdydIdy
  ConvolveX(dIdydIdy, mask, gaussI, vvftmp);
  ConvolveY(dIdydIdy, mask, gaussI, vvftmp);

  //----------------------------------------------------------------------
  // dIdxdIdy
  ConvolveX(dIdxdIdy, mask, gaussI, vvftmp);
  ConvolveY(dIdxdIdy, mask, gaussI, vvftmp);
}

void Harris::preprocess(void) {
  vector<vector<Vec3f>> vvvftmp;
  vvvftmp.resize(height);
  for (int y = 0; y < height; ++y) {
    vvvftmp[y].resize(width);
    for (int x = 0; x < width; ++x)
      vvvftmp[y][x] = Vec3f();
  }

  dIdx = image;

  vector<float> dfilter, ifilter;
  dfilter.resize(3);
  dfilter[0] = -0.5;
  dfilter[1] = 0;
  dfilter[2] = 0.5;
  ifilter.resize(3);
  ifilter[0] = 1.0 / 3.0;
  ifilter[1] = 1.0 / 3.0;
  ifilter[2] = 1.0 / 3.0;

  ConvolveX(dIdx, mask, dfilter, vvvftmp);
  ConvolveY(dIdx, mask, ifilter, vvvftmp);

  dIdy = image;
  ConvolveX(dIdy, mask, ifilter, vvvftmp);
  ConvolveY(dIdy, mask, dfilter, vvvftmp);
}

void Harris::setResponse(void) {
  response.clear();
  response.resize(height);
  for (int y = 0; y < height; ++y) {
    response[y].resize(width);
    for (int x = 0; x < width; ++x) {
      response[y][x] = 0.0;
      if (!mask.empty() && !mask[y][x])
        continue;

      const float D  = dIdxdIdx[y][x] * dIdydIdy[y][x] -
                       dIdxdIdy[y][x] * dIdxdIdy[y][x];
      const float tr = dIdxdIdx[y][x] + dIdydIdy[y][x];
      response[y][x] = D - 0.06 * tr * tr;
    }
  }

  //----------------------------------------------------------------------
  // suppress non local max
  vector<vector<float>> vvftmp = response;
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 1; x < width - 1; ++x) {
      if (response[y][x] < response[y    ][x + 1] ||
          response[y][x] < response[y    ][x - 1] ||
          response[y][x] < response[y + 1][x    ] ||
          response[y][x] < response[y - 1][x    ])
        vvftmp[y][x] = 0.0;
    }
  }

  vvftmp.swap(response);
}

void Harris::Run(const std::vector<unsigned char> &image,
                 const std::vector<unsigned char> &mask,
                 const std::vector<unsigned char> &edge, const int width,
                 const int height, const int gspeedup, const float sigma,
                 std::multiset<Point> &result) {

  cerr << "Harris running ..." << flush;
  this->width  = width;
  this->height = height;
  sigmaD = sigma;
  sigmaI = sigma;
  init(image, mask, edge);
  setDerivatives();
  setResponse();

  const int factor = 2;
  const int maxPointsGrid = factor * factor;
  const int gridsize = gspeedup * factor;

  const int w = (width  + gridsize - 1) / gridsize;
  const int h = (height + gridsize - 1) / gridsize;

  vector<vector<multiset<Point>>> resultgrids;
  resultgrids.resize(h);
  for (int y = 0; y < h; ++y)
    resultgrids[y].resize(w);

  const int margin = (int)gaussD.size() / 2;
  for (int y = margin; y < height - margin; ++y) {
    for (int x = margin; x < width - margin; ++x) {
      if (response[y][x] == 0.0)
        continue;

      const int x0 = min(x / gridsize, w - 1);
      const int y0 = min(y / gridsize, h - 1);

      if ((int)resultgrids[y0][x0].size() < maxPointsGrid || resultgrids[y0][x0].begin()->response < response[y][x]) {
        Point p;
        p.iCoord = Vec3f(x, y, 1.0f);
        p.response = response[y][x];
        p.type = 0;

        resultgrids[y0][x0].insert(p);
        if (maxPointsGrid < (int)resultgrids[y0][x0].size())
          resultgrids[y0][x0].erase(resultgrids[y0][x0].begin());
      }
    }
  }

  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x) {
      // const float threshold = setThreshold(resultgrids[y][x]);
      multiset<Point>::iterator begin = resultgrids[y][x].begin();
      multiset<Point>::iterator end   = resultgrids[y][x].end();
      while (begin != end) {
        // if (threshold <= begin->response)
        result.insert(*begin);
        begin++;
      }
    }

  cerr << (int)result.size() << " harris done" << endl;
}
