// MedianVariance.h

#ifndef MEDIANVARIANCE_H

#define MEDIANVARIANCE_H

#include <algorithm>
#include <vector>
#include "Variance.h"

template <class T> class MedianVariance {
public:
  static constexpr int PERCENTILE = 25;
public:
  MedianVariance(T first=0, int chunksize0=250):
    col(first), chunksize(chunksize0) { }
  void reset(T first=0) {
    col.reset(first);
    vars.erase(vars.begin(),vars.end());
    means.erase(means.begin(), means.end());
    i = chunksize;
  }
  void addexample(T d) {
    col.addexample(d);
    if (!--i) {
      means.push_back(col.mean());
      vars.push_back(col.var());
      col.reset(col.mean());
      i=chunksize;
    }
  }
  T mean() { // actually: median of means
    int n=means.size()*50/100;
    if (n==0) return 0;
    typename std::vector<T>::iterator i=means.begin()+n;
    std::nth_element(means.begin(),i,means.end());
    return *i;
  }
  T var() { // actually: PERCENTILE of vars
    int n=vars.size()*PERCENTILE/100;
    if (n==0)
      return 0;
    auto i = vars.begin()+n;
    std::nth_element(vars.begin(),i,vars.end());
    return *i;
  }
  int chunks() const { return vars.size(); }
private:
  Variance<T> col;
  std::vector<T> means;
  std::vector<T> vars;
  int chunksize;
  int i;
};



#endif
