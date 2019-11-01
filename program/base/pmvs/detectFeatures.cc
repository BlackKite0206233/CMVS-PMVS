#include "detectFeatures.h"
#include "../image/image.h"
#include "dog.h"
#include "harris.h"
#include "point.h"
#include <fstream>
#include <iostream>

using namespace PMVS3;
using namespace std;
using namespace img;

DetectFeatures::DetectFeatures(void) {
  mtx_init(&rwLock, mtx_plain | mtx_recursive);
}

DetectFeatures::~DetectFeatures() { mtx_destroy(&rwLock); }

void DetectFeatures::Run(const PhotoSet &pss, const int num, const int csize, const int level, const int CPU) {
  ps          = &pss;
  cSize       = csize;
  this->level = level;
  this->CPU   = CPU;

  points.clear();
  points.resize(num);

  //----------------------------------------------------------------------
  for (int index = 0; index < num; ++index)
    jobs.push_back(index);

  vector<thrd_t> threads(this->CPU);
  for (auto& t : threads)
    thrd_create(&t, &runThreadTmp, (void *)this);
  for (auto& t : threads)
    thrd_join(t, NULL);
  //----------------------------------------------------------------------
  cerr << "done" << endl;
}

int DetectFeatures::runThreadTmp(void *arg) {
  DetectFeatures *detectFeatures = (DetectFeatures *)arg;
  detectFeatures->runThread();
  return 0;
}

void DetectFeatures::runThread(void) {
  while (1) {
    int index = -1;
    mtx_lock(&rwLock);
    if (!jobs.empty()) {
      index = jobs.front();
      jobs.pop_front();
    }
    mtx_unlock(&rwLock);
    if (index == -1)
      break;

    const int image = ps->images[index];
    cerr << image << ' ' << flush;

    //?????????????  May need file lock, because targetting images
    // should not overlap among multiple processors.
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.affin%d", ps->prefix.c_str(), image, level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (ifstr.is_open()) {
      ifstr.close();
      continue;
    }
    ifstr.close();

    //----------------------------------------------------------------------
    // parameters
    // for harris
    const float sigma = 4.0f;
    // for DoG
    const float firstScale = 1.0f;
    const float lastScale  = 3.0f;

    //----------------------------------------------------------------------
    // Harris
    {
      Harris harris;
      multiset<Point> result;
      harris.Run(ps->photos[index].GetImage(level),
                 ps->photos[index].Image::GetMask(level),
                 ps->photos[index].Image::GetEdge(level),
                 ps->photos[index].GetWidth(level),
                 ps->photos[index].GetHeight(level), cSize, sigma,
                 result);

			for (auto& r : result) 
				points[index].push_back(r);
    }

    //----------------------------------------------------------------------
    // DoG
    {
      Dog dog;
      multiset<Point> result;
      dog.Run(ps->photos[index].GetImage(level),
              ps->photos[index].Image::GetMask(level),
              ps->photos[index].Image::GetEdge(level),
              ps->photos[index].GetWidth(level),
              ps->photos[index].GetHeight(level), cSize, firstScale,
              lastScale, result);

			for (auto& r : result) 
				points[index].push_back(r);
    }
  }
}
