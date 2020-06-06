// NoiseLevels.h

#ifndef NOISELEVELS_H

#define NOISELEVELS_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "MedianVariance.h"

class NoiseLevels {
public:
  static constexpr int MINCHUNKS = 5;
  NoiseLevels() {
    reset();
  }
  double mean() const {
    return means;
  }
  double std() const {
    return stds;
  }
  void reset() {
    ready = false;
    mv.reset();
  }
  void train(CyclBuf<raw_t> const &buf, int start, int end) {
    for (int k=start; k<end; k++)
      mv.addexample(buf[k]);
  }
  void crash(char const *msg) {
    std::cerr << msg << "\n";
    std::exit(1);
  }
  void makeready() {
    if (chunks() < MINCHUNKS)
      crash("Too few chunks to compute meaningful estimates");
    means = mv.mean();
    stds = std::sqrt(mv.var());
    ready = true;
  }
  int chunks() const {
    return mv.chunks();
  }
  bool isready() const {
    return ready;
  }
private:
  MedianVariance<double> mv;
  bool ready;
  double means;
  double stds;
};

#endif
