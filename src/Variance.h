// Variance.h

#ifndef VARIANCE_H

#define VARIANCE_H

template <class T> class Variance {
public:
  Variance(T first=0) { reset(first); }
  void reset(T first=0) { av0=first; sx=0; sxx=0; n=0; }
  void addexample(T d) { d-=av0; sx+=d; sxx+=d*d; n++; }
  T mean() const { return av0 + sx/n; }
  T var() const { return (sxx-sx*sx/n)/(n-1); }
private:
  T av0, sx, sxx;
  int n;
};

#endif
