#include "bmap.h"
#include <cassert>

BMap::BMap(RangeCoderSH &rc)
:rc(rc),pmix(4),
ctx0(0),ctx1(0),ctx2(0),ctx3(0),
pv(4)
{
}

int32_t BMap::MirrorCtx(const BIntMap &bmap,int32_t idx) const
{
  int ctx=0;

  int64_t v=bmap.idx2val(idx); //index to val
  int64_t mv=-v;//mirror at origin
  int64_t mi=bmap.val2idx(mv);//val to index
  static constexpr int ofs[]={0,-1,1,-2,2};
  for (std::size_t k=0;k<std::size(ofs);k++) {
    int64_t j=mi+ofs[k];
    int b=0;
    if (j>=0 && j<(int)bmap.buf.size() && j<idx)
      b=bmap.buf[j]?1:0;
    ctx|=(b<<k);
  }
  return ctx;
}

int32_t BMap::RefCtx(const BIntMap &bmap_ref,int32_t val) const
{
  std::size_t bsize=bmap_ref.buf.size();
  int32_t ctx=0;
  if (bsize) {
    int mi=bmap_ref.val2idx(val);
    static constexpr int ofs[]={0,-1,1,-2,2};
    for (std::size_t k=0;k<std::size(ofs);k++) {
      int64_t j=mi+ofs[k];
      int b=0;
      if (j>=0 && j<(int)bsize) {
        b=bmap_ref.buf[j]?1:0;
      }
      ctx|=(b<<k);
    }
  }
  return ctx;
}

void BMap::Update(int bit)
{
  c0[ctx0].update(bit,WCNT);
  c1[ctx1].update(bit,WCNT);
  c2[ctx2].update(bit,WCNT);
  c3[ctx3].update(bit,WCNT);
  sse.Update(bit,WSSE);
  pmix.Update(bit,WMIX);
}

int32_t BMap::Predict()
{
  pv[0]=c0[ctx0].p1;
  pv[1]=c1[ctx1].p1;
  pv[2]=c2[ctx2].p1;
  pv[3]=c3[ctx3].p1;

  int32_t pm=pmix.Predict(pv);
  int32_t p_sse = sse.Predict(pm);
  int32_t px=(3*pm+p_sse)>>2;
  return px;
}

void BMap::Encode(const BIntMap &bmap,const BIntMap &)
{
  std::size_t buf_len = bmap.buf.size();
  assert(buf_len==(bmap.maxval-bmap.minval+1));
  //std::cout << "bmap (" << buf_len << "): " << bmap.minval << ' ' << bmap.maxval << '\n';

  //linear scan order
  RunStats rs;
  for (int i=0;i<(int)buf_len;i++)
  {
    ctx0 = rs.bhist; //history
    ctx1 = (rs.rt<<1)|rs.rb; //run
    ctx2 = MirrorCtx(bmap,i); //mirror
    ctx3 = 0;//RefCtx(bmap_ref,bmap.idx2val(i)); //reference

    int32_t px=Predict();
    const int bit=bmap.buf[i]?1:0;
    rc.EncodeBitOne(px,bit);

    Update(bit);

    rs.Update(bit);
  }
}

void BMap::Decode(BIntMap &bmap,const BIntMap &)
{
  std::size_t buf_len = bmap.buf.size();
  if (buf_len!=static_cast<std::size_t>(bmap.maxval-bmap.minval+1))
    throw std::runtime_error("BMap: wrong allocation");

  //std::cout << "bmap (" << buf_len << "): " << bmap.minval << ' ' << bmap.maxval << '\n';

  //linear scan order
  RunStats rs;
  for (int i=0;i<(int)buf_len;i++)
  {
    ctx0 = rs.bhist; //history
    ctx1 = (rs.rt<<1)|rs.rb; //run
    ctx2 = MirrorCtx(bmap,i); //mirror
    ctx3 = 0;//RefCtx(bmap_ref,bmap.idx2val(i)); //reference

    int32_t px=Predict();
    int bit=rc.DecodeBitOne(px);

    Update(bit);
    rs.Update(bit);

    bmap.buf[i] = bit;
  }
}

