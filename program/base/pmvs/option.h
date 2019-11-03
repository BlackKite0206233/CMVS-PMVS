#ifndef PMVS3_OPTION_H
#define PMVS3_OPTION_H

#include <map>
#include <string>
#include <vector>

namespace PMVS3 {

struct Option {
public:
  int   level;
  int   cSize;
  float threshold;
  int   wSize;
  int   minImageNum;
  int   CPU;
  float setEdge;
  int   useBound;
  int   useVisData;
  int   sequence;

  float maxAngleThreshold;
  float quadThreshold;

  std::string prefix;
  std::string option;
  int tFlag;
  std::vector<int> tImages;
  int oFlag;
  std::vector<int> oImages;

  std::map<int, int> dict;

  std::vector<int> bindexes;
  std::vector<std::vector<int>> visData;
  std::vector<std::vector<int>> visData2;

  Option(void);

  void Init(const std::string prefix, const std::string option);

protected:
  void initOimages(void);
  void initBindexes(const std::string sbimages);
  void initVisdata(void);
  void initVisdata2(void);
};
}; // namespace PMVS3

#endif // PMVS3_OPTION_H
