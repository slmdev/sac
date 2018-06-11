#ifndef HISTBUF_H
#define HISTBUF_H

#include "../global.h"

// rolling buffer, n must be power of two
template <class T>
class HistBuffer {
  public:
      HistBuffer(int n)
      :n(n),mask(n-1),pos(0),buf(n)
      {
      }
      T& operator[](std::size_t idx) { return buf[(pos-idx)&mask];};
      const T operator[](std::size_t idx) const { return buf[(pos-idx)&mask];};
      void PushBack(T val)
      {
         buf[pos]=val;
         pos=(pos+1)&mask;
      }
  private:
    int n,mask,pos;
    std::vector <T> buf;
};


#endif // HISTBUF_H
