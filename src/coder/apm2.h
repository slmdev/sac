#ifndef APM2_H
#define APM2_H

//prediction as correction to p1
//using LogMixerBias and multiple sse-layer
class APM2 {
  static constexpr float DCONF = 0.15;
  public:
    //nmix=mixer context size
    //nctx=vector of sse table sizes
    //mconf=initial mixer confidence (0.0=pass-through, 1.0=uniform
    //fscale=saturation scale for sse, scale=myDomain.Fwd(PSCALE) for fscale=1.0
    APM2(int nmix,const std::vector<std::size_t>&nctx,float mconf=0.0,float fscale=1.0)
    :n(nctx.size()),
     sse_tbl(n),
     sse_mix(nmix,LogMixerBias(n+2)),
     input(n+2),
     psse(n)
    {
      for (std::size_t i=0;i<n;i++) {
        //nctx[i] constructed SSE tables
        sse_tbl[i].resize(
          nctx[i],SSENL<16>(fscale));
      }
      if (mconf>0.0) {
        for (auto &m:sse_mix)
          m.SetConfidence(mconf);
      }
    }
    int Predict(int p1,int mix_ctx,std::span<const int>ctx)
    {
      int64_t sum=0;
      for (std::size_t i=0;i<n;i++) {
        psse[i]=&sse_tbl[i][ctx[i]];
        //every sse table gets same p1
        input[i] = psse[i]->Predict(p1);
        sum+=input[i];
      }
      //avg sse prediction
      const int p_avg = static_cast<int>((sum+(n/2))/n);
      const int l_avg = myDomain.Fwd(p_avg);
      //less confidence
      input[n] = myDomain.Inv((1.0-DCONF)*l_avg);
      //more confidence
      input[n+1] = myDomain.Inv((1.0+DCONF)*l_avg);
      pmix = &sse_mix[mix_ctx];
      return pmix->Predict(p1,input);
    }
    void Update(int bit,int climit,int mrate)
    {
      for (auto* tbl:psse)
        tbl->Update(bit,climit);

      pmix->Update(bit,mrate);
    }
  private:
    std::size_t n;
    //2D array
    std::vector<std::vector<SSENL<16>>>sse_tbl;
    std::vector<LogMixerBias> sse_mix;

    std::vector<int> input;
    std::vector<SSENL<16>*> psse;
    LogMixerBias *pmix=nullptr;
};

#endif

