// CyclBuf.h

#ifndef CYCLBUF_H

#define CYCLBUF_H

#include <vector>
#include <cstdint>

template <class T> class CyclBuf {
 public:
  CyclBuf(int log2size=16):
    stride(1) {
    int size = 1 << log2size;
    mask = size - 1;
    vec = std::vector<T>(size, 0);
    data = vec.data();
  }
  CyclBuf(T *data, int log2size, int stride=1):
    data(data), stride(stride) {
    int size = 1 << log2size;
    mask = size - 1;
  }
  T const &operator[](std::uint32_t index) const {
    index &= mask;
    return data[index*stride];
  }
  T &operator[](std::uint32_t index) {
    index &= mask;
    return data[index*stride];
  }
 private:
  std::vector<T> vec;
  T *data;
  std::uint32_t stride;
  std::uint32_t mask;
};

#endif
