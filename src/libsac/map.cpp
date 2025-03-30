#include "map.h"

MapEncoder::MapEncoder(RangeCoderSH &rc,std::vector <bool>&usedl,std::vector <bool>&usedh)
:rc(rc),mixl(4,NMixLogistic(5)),mixh(4,NMixLogistic(5)),finalmix(2),ul(usedl),uh(usedh)
{
}

int MapEncoder::PredictLow(int i)
{
  int ctx1=ul[i-1];
  int ctx2=uh[i-1];
  int ctx3=i>1?ul[i-2]:0;

  pc1=&cnt[ctx1];
  pc2=&cnt[2+ctx2];
  pc3=&cnt[4+(ctx1<<1)+ctx3];
  pc4=&cnt[8+(ctx1<<1)+ctx2];

  int sctx=ul[i-1];
  if (i>1) sctx+=(ul[i-2]<<1);
  if (i>2) sctx+=(ul[i-3]<<2);
  if (i>3) sctx+=(ul[i-4]<<3);
  px=&cctx[sctx];

  mix=&mixl[ctx1+(ctx3<<1)];
  std::vector <int>p={pc1->p1,pc2->p1,pc3->p1,pc4->p1,px->p1};
  return mix->Predict(p);
}

int MapEncoder::PredictHigh(int i)
{
  int ctx1=uh[i-1];
  int ctx2=ul[i];
  int ctx3=i>1?uh[i-2]:0;
  //int n=0;
  pc1=&cnt[12+ctx1];
  pc2=&cnt[12+2+ctx2];
  pc3=&cnt[12+4+(ctx1<<1)+ctx3];
  pc4=&cnt[12+8+(ctx1<<1)+ctx2];

  int sctx=uh[i-1];
  if (i>1) sctx+=(uh[i-2]<<1);
  if (i>2) sctx+=(uh[i-3]<<2);
  if (i>3) sctx+=(uh[i-4]<<3);
  px=&cctx[32+sctx];
  mix=&mixh[ctx1+(ctx3<<1)];
  std::vector <int>p={pc1->p1,pc2->p1,pc3->p1,pc4->p1,px->p1};
  return mix->Predict(p);
}

void MapEncoder::Update(int bit)
{
  pc1->update(bit,cnt_upd_rate);
  pc2->update(bit,cnt_upd_rate);
  pc3->update(bit,cnt_upd_rate);
  pc4->update(bit,cnt_upd_rate);
  px->update(bit,cnt_upd_rate);
  mix->Update(bit,mix_upd_rate);
}

int MapEncoder::PredictSSE(int p1,int ctx)
{
  std::vector <int>vp={sse[ctx].Predict(p1),p1};
  return finalmix.Predict(vp);
}

void MapEncoder::UpdateSSE(int bit,int ctx)
{
  sse[ctx].Update(bit,cntsse_upd_rate);
  finalmix.Update(bit,mixsse_upd_rate);
}

void MapEncoder::Encode()
{
  for (int i=1;i<=1<<15;i++) {
    int bit=ul[i];

    rc.EncodeBitOne(PredictSSE(PredictLow(i),0),bit);
    Update(bit);
    UpdateSSE(bit,0);

    bit=uh[i];
    rc.EncodeBitOne(PredictSSE(PredictHigh(i),0),bit);
    Update(bit);
    UpdateSSE(bit,0);
  }
}

void MapEncoder::Decode()
{
  for (int i=1;i<=1<<15;i++) {
    int bit=rc.DecodeBitOne(PredictSSE(PredictLow(i),0));
    Update(bit);
    ul[i]=bit;
    UpdateSSE(bit,0);

    bit=rc.DecodeBitOne(PredictSSE(PredictHigh(i),0));
    Update(bit);
    uh[i]=bit;
    UpdateSSE(bit,0);
  }
}

Remap::Remap()
:scale(1<<15),usedl(scale+1),usedh(scale+1)
{
}

void Remap::Reset()
{
  std::fill(begin(usedl),end(usedl),0);
  std::fill(begin(usedh),end(usedh),0);
  vmin=vmax=0;
}

double Remap::Compare(const Remap &cmap)
{
  int diff=0;
  for (int i=1;i<=scale;i++) {
    if (usedl[i]!=cmap.usedl[i]) diff++;
    if (usedh[i]!=cmap.usedh[i]) diff++;
  }
  return diff*100./double(2*scale);
}

void Remap::Analyse(int32_t *src,int numsamples)
{
  for (int i=0;i<numsamples;i++) {
    int val=src[i];
    if (val>0) {
      if (val>scale) std::cout << "val too large: " << val << std::endl;
      else {
        if (val>vmax) vmax=val;
        usedh[val]=true;
      }
    } else if (val<0) {
      val=(-val);
      if (val>scale) std::cout << "val too large: " << val << std::endl;
      else {
        if (val>vmin) vmin=val;
        usedl[val]=true;
      }
    }
  }
  mapl.resize((1<<15)+1);
  maph.resize((1<<15)+1);
  int j=1;
  for (int i=1;i<=(1<<15);i++) {
    mapl[i]=j;
    if (usedl[i]) {j++;};
  }
  j=1;
  for (int i=1;i<=(1<<15);i++) {
    maph[i]=j;
    if (usedh[i]) {j++;};
  }
}

bool Remap::isUsed(int val)
{
  if (val>scale) return false;
  if (val<-scale) return false;
  if (val>0) return usedh[val];
  if (val<0) return usedl[-val];
  return true;
}

int32_t Remap::Map2(int32_t pred)
{
  if (pred>0) return maph[pred];
  else if (pred<0) return -mapl[-pred];
  else return 0;
}

int32_t Remap::Map(int32_t pred,int32_t err)
{
  int sgn=1;
  if (err==0) return 0;
  if (err<0) {err=-err;sgn=-1;};

  int merr=0;
  for (int i=1;i<=err;i++) {
    if (isUsed(pred+(sgn*i))) merr++;
  }
  return sgn*merr;
}

int32_t Remap::Unmap(int32_t pred,int32_t merr)
{
  int sgn=1;
  if (merr==0) return 0;
  if (merr<0) {merr=-merr;sgn=-1;};

  int err=1;
  int terr=0;
  while (1) {
    if (isUsed(pred+(sgn*err))) terr++;
    if (terr==merr) break;
    err++;
  }
  return sgn*err;
}
