#include "../base/pmvs/findMatch.h"
#include "../base/pmvs/option.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace PMVS3;
using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " prefix option_file [Optional export]" << endl
         << endl
         << "--------------------------------------------------" << endl
         << "level       1    csize    2" << endl
         << "threshold   0.7  wsize    7" << endl
         << "minImageNum 3    CPU      4" << endl
         << "useVisData  0    sequence -1" << endl
         << "quad        2.5  maxAngle 10.0" << endl
         << "--------------------------------------------------" << endl
         << "2 ways to specify targetting images" << endl
         << "timages  5  1 3 5 7 9 (enumeration)" << endl
         << "        -1  0 24 (range specification)" << endl
         << "--------------------------------------------------" << endl
         << "4 ways to specify other images" << endl
         << "oimages  5  0 2 4 6 8 (enumeration)" << endl
         << "        -1  24 48 (range specification)" << endl
         << endl
         << "[Optional export] PATCH PSET" << endl
         << " i.e export patch and pset: prefix option_file PATCH PSET"
         << " i.e export patch only: prefix option_file PATCH" << endl;
    exit(1);
  }

  for (int i = 0; i < argc; ++i) {
    cout << endl << argv[i];
  }
  cout << std::endl;

  clock_t begin = clock();

  PMVS3::Option    option;
  PMVS3::FindMatch findMatch;
  {
	option.Init(argv[1], argv[2]);

	findMatch.Init(option);
	findMatch.Run();
  }

  cout << "total time: " << (double)(clock() - begin) / CLOCKS_PER_SEC << endl;

  bool bExportPLY   = true;
  bool bExportPatch = false;
  bool bExportPSet  = false;

  for (int i = 3; i < argc; ++i) {
    std::string option(argv[i]);
    if (option == "PATCH")
      bExportPatch = true;
    if (option == "PSET")
      bExportPSet = true;
  }

	stringstream ss;
	string str;
	ss << argv[1] << "models/" << argv[2];
	ss >> str;
  findMatch.Write(str, bExportPLY, bExportPatch, bExportPSet);
}
