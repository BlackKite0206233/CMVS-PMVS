#include <cstdlib>
#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include "option.h"
#include <algorithm>
#include <math.h>

using namespace std;
using namespace PMVS3;

Option::Option(void) {
  level       = 1;
  cSize       = 2;
  threshold   = 0.7;
  wSize       = 7;
  minImageNum = 3;
  CPU         = 4;
  setEdge     = 0.0f;
  useBound    = 0;
  useVisData  = 0;
  sequence    = -1;
  tFlag       = -10;
  oFlag       = -10;

  // Max angle must be at least this big
  maxAngleThreshold = 10.0f * M_PI / 180.0f;
  // The smaller the tighter
  quadThreshold = 2.5f;
}

void Option::Init(const std::string prefix, const std::string option) {
  this->prefix = prefix;
  this->option = option;
  std::ifstream ifstr;
  string optionfile = this->prefix + this->option;
  ifstr.open(optionfile.c_str());
  while (1) {
    string name;
    ifstr >> name;
    if (ifstr.eof())
      break;
    if (name[0] == '#') {
      char buffer[1024];
      ifstr.putback('#');
      ifstr.getline(buffer, 1024);
      continue;
    }
    if (name == "level")
      ifstr >> level;
    else if (name == "csize")
      ifstr >> cSize;
    else if (name == "threshold")
      ifstr >> threshold;
    else if (name == "wsize")
      ifstr >> wSize;
    else if (name == "minImageNum")
      ifstr >> minImageNum;
    else if (name == "CPU")
      ifstr >> CPU;
    else if (name == "setEdge")
      ifstr >> setEdge;
    else if (name == "useBound")
      ifstr >> useBound;
    else if (name == "useVisData")
      ifstr >> useVisData;
    else if (name == "sequence")
      ifstr >> sequence;
    else if (name == "timages") {
      ifstr >> tFlag;
      if (tFlag == -1) {
        int firstimage, lastimage;
        ifstr >> firstimage >> lastimage;
        for (int i = firstimage; i < lastimage; ++i)
          tImages.push_back(i);
      } else if (0 < tFlag) {
        for (int i = 0; i < tFlag; ++i) {
          int itmp;
          ifstr >> itmp;
          tImages.push_back(itmp);
        }
      } else {
        cerr << "tflag is not valid: " << tFlag << endl;
        exit(1);
      }
    } else if (name == "oimages") {
      ifstr >> oFlag;
      if (oFlag == -1) {
        int firstimage, lastimage;
        ifstr >> firstimage >> lastimage;
        for (int i = firstimage; i < lastimage; ++i)
          oImages.push_back(i);
      } else if (0 <= oFlag) {
        for (int i = 0; i < oFlag; ++i) {
          int itmp;
          ifstr >> itmp;
          oImages.push_back(itmp);
        }
      } else if (oFlag != -2 && oFlag != -3) {
        cerr << "oflag is not valid: " << oFlag << endl;
        exit(1);
      }
    } else if (name == "quad")
      ifstr >> quadThreshold;
    else if (name == "maxAngle") {
      ifstr >> maxAngleThreshold;
      maxAngleThreshold *= M_PI / 180.0f;
    } else {
      cerr << "Unrecognizable option: " << name << endl;
      exit(1);
    }
  }
  ifstr.close();

  if (tFlag == -10 || oFlag == -10) {
    cerr << "tflag and oflag not specified: " << tFlag << ' ' << oFlag << endl;
    exit(1);
  }

  //----------------------------------------------------------------------
  string sbimages = prefix + string("bimages.dat");

  for (int i = 0; i < (int)tImages.size(); ++i)
    dict[tImages[i]] = i;

  initOimages();
  initVisdata();

  if (useBound)
    initBindexes(sbimages);

  cerr << "--------------------------------------------------" << endl
       << "--- Summary of specified options ---" << endl;
  cerr << "# of timages: " << (int)tImages.size();
  if (tFlag == -1)
    cerr << " (range specification)" << endl;
  else
    cerr << " (enumeration)" << endl;
  cerr << "# of oimages: " << (int)oImages.size();
  if (oFlag == -1)
    cerr << " (range specification)" << endl;
  else if (0 <= oFlag)
    cerr << " (enumeration)" << endl;
  else if (oFlag == -2)
    cerr << " (vis.dat is used)" << endl;
  else if (oFlag == -3)
    cerr << " (not used)" << endl;

  cerr << "level: " << level << "  csize: " << cSize << endl
       << "threshold: " << threshold << "  wsize: " << wSize << endl
       << "minImageNum: " << minImageNum << "  CPU: " << CPU << endl
       << "useVisData: " << useVisData << "  sequence: " << sequence << endl;
  cerr << "--------------------------------------------------" << endl;
}

void Option::initOimages(void) {
  if (oFlag != -2)
    return;

  string svisData = prefix + string("vis.dat");
  ifstream ifstr;
  ifstr.open(svisData.c_str());
  if (!ifstr.is_open()) {
    cerr << "No vis.dat although specified to initOimages: " << endl
         << svisData << endl;
    exit(1);
  }

  string header;
  int num2;
  ifstr >> header >> num2;

  oImages.clear();
  for (int c = 0; c < num2; ++c) {
    int index0;
    map<int, int>::iterator ite0 = dict.find(c);
    if (ite0 == dict.end())
      index0 = -1;
    else
      index0 = ite0->second;
    int itmp;
    ifstr >> itmp >> itmp;
    for (int i = 0; i < itmp; ++i) {
      int itmp2;
      ifstr >> itmp2;
      if (index0 != -1 && dict.find(itmp2) == dict.end())
        oImages.push_back(itmp2);
    }
  }
  ifstr.close();

  sort(oImages.begin(), oImages.end());
  oImages.erase(unique(oImages.begin(), oImages.end()), oImages.end());
}

// When do not use vis.dat
void Option::initVisdata(void) {
  // Case classifications. Set visData by using vis.dat or not.
  if (useVisData == 0) {
    const int tnum = (int)tImages.size();
    const int onum = (int)oImages.size();
    const int num  = tnum + onum;
    visData.resize(num);
    visData2.resize(num);
    for (int y = 0; y < num; ++y) {
      visData[y].resize(num);
      for (int x = 0; x < num; ++x)
        if (x == y)
          visData[y][x] = 0;
        else {
          visData[y][x] = 1;
          visData2[y].push_back(x);
        }
    }
  } else
    initVisdata2();
}

// Given timages and oimages, set visData, visData2
void Option::initVisdata2(void) {
  string svisData = prefix + string("vis.dat");

  vector<int> images;
  images.insert(images.end(), tImages.begin(), tImages.end());
  images.insert(images.end(), oImages.begin(), oImages.end());
  map<int, int> dict2;
  for (int i = 0; i < (int)images.size(); ++i)
    dict2[images[i]] = i;

  ifstream ifstr;
  ifstr.open(svisData.c_str());
  if (!ifstr.is_open()) {
    cerr << "No vis.dat although specified to initVisdata2: " << endl
         << svisData << endl;
    exit(1);
  }

  string header;
  int num2;
  ifstr >> header >> num2;

  visData2.resize((int)images.size());
  for (int c = 0; c < num2; ++c) {
    int index0;
    map<int, int>::iterator ite0 = dict2.find(c);
    if (ite0 == dict2.end())
      index0 = -1;
    else
      index0 = ite0->second;
    int itmp;
    ifstr >> itmp >> itmp;
    for (int i = 0; i < itmp; ++i) {
      int itmp2;
      ifstr >> itmp2;
      int index1;
      map<int, int>::iterator ite1 = dict2.find(itmp2);
      if (ite1 == dict2.end())
        index1 = -1;
      else
        index1 = ite1->second;

      if (index0 != -1 && index1 != -1)
        visData2[index0].push_back(index1);
    }
  }
  ifstr.close();

  const int num = (int)images.size();
  visData.clear();
  visData.resize(num);
  for (int y = 0; y < num; ++y) {
    visData[y].resize(num);
    fill(visData[y].begin(), visData[y].end(), 0);
    for (int x = 0; x < (int)visData2[y].size(); ++x)
      visData[y][visData2[y][x]] = 1;
  }

  // check symmetry
  for (int i = 0; i < (int)visData.size(); ++i) {
    for (int j = i + 1; j < (int)visData.size(); ++j) {
      if (visData[i][j] != visData[j][i]) {
        visData[i][j] = visData[j][i] = 1;
      }
    }
  }
}

void Option::initBindexes(const std::string sbimages) {
  if (sbimages.empty())
    return;

  bindexes.clear();
  ifstream ifstr;
  ifstr.open(sbimages.c_str());
  if (!ifstr.is_open()) {
    cerr << "File not found: " << sbimages << endl;
    exit(1);
  }

  cerr << "Reading bimages" << endl;
  int itmp;
  ifstr >> itmp;
  for (int i = 0; i < itmp; ++i) {
    int itmp0;
    ifstr >> itmp0;

    if (dict.find(itmp0) != dict.end())
      bindexes.push_back(dict[itmp0]);
  }
  ifstr.close();
}
