#ifndef BUFIO_H
#define BUFIO_H

#include <vector>
#include <cstdint>

class BufIO {
  public:
      BufIO():buf(1024){Reset();};
      BufIO(int initsize):buf(initsize){Reset();};
      void Reset(){bufpos=0;};
      void PutByte(int val)
      {
        if (bufpos>=buf.size()) buf.resize(buf.size()*2);
        buf[bufpos++]=val;
      }
      int GetByte() {
        if (bufpos>=buf.size()) return -1;
        else return buf[bufpos++];
      }
      size_t GetBufPos(){return bufpos;};
      std::vector <uint8_t> &GetBuf(){return buf;};
  private:
     size_t bufpos;
     std::vector <uint8_t>buf;
};

#endif
