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
      T& operator[](int idx) { return buf[(pos-idx)&mask];};
      const T operator[](int idx) const { return buf[(pos-idx)&mask];};
      void PushBack(T val)
      {
         buf[pos]=val;
         pos=(pos+1)&mask;
      }
  private:
    int n,mask,pos;
    std::vector <T> buf;
};

template <class T>
class RollBuffer {
  public:
    RollBuffer(int n)
    :n(n),pos(-1),buf(n)
    {

    }
    const T operator[](int idx) const
    {
      return buf[clamp_idx(pos-idx)];
    };
    void PushBack(T val)
    {
      if (++pos>=n) pos=0;
      buf[pos]=val;
    }
    int clamp_idx(int idx) const
    {
      if (idx>=n) idx-=n;
      else if (idx<0) idx+=n;
      return idx;
    }
    const std::vector<T> &getbuf() {return buf;};
  private:
    int n,pos;
    std::vector <T> buf;
};

#endif // HISTBUF_H
