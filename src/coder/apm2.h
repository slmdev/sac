#ifndef APM2_H
#define APM2_H

//prediction as correction to p1
//using LogMixerBias and multiple sse-layer
class APM2 {
  static constexpr float DCONF = 0.15;
  public:
    APM2(int nmix,const std::vector<std::size_t>&nctx)
    :n(nctx.size()),
     sse_tbl(n),
     sse_mix(nmix,LogMixerBias(n+2)),
     input(n+2),
     psse_tbl(n)
    {
      for (std::size_t i=0;i<n;i++) {
        //nctx[i] default-constructed SSE tables
        sse_tbl[i].resize(nctx[i]);
      }
    }
    int Predict(int p1,int mix_ctx,std::span<const int>ctx)
    {
      pmix_ctx = &sse_mix[mix_ctx];
      int64_t sum=0;
      for (std::size_t i=0;i<n;i++) {
        psse_tbl[i]=&sse_tbl[i][ctx[i]];
        //every sse table gets same p1
        input[i] = psse_tbl[i]->Predict(p1);
        sum+=input[i];
      }
      //avg sse prediction
      const int p_avg = static_cast<int>((sum+(n/2))/n);
      const int l_avg = myDomain.Fwd(p_avg);
      //less confidence
      input[n] = myDomain.Inv((1.0-DCONF)*l_avg);
      //more confidence
      input[n+1] = myDomain.Inv((1.0+DCONF)*l_avg);
      return sse_mix[mix_ctx].Predict(p1,input);
    }
    void Update(int bit,int climit,int mrate)
    {
      for (auto* tbl:psse_tbl)
        tbl->Update(bit,climit);

      pmix_ctx->Update(bit,mrate);
    }
  private:
    std::size_t n;
    //2D array
    std::vector<std::vector<SSENL<16>>>sse_tbl;
    std::vector<LogMixerBias> sse_mix;

    std::vector<int> input;
    std::vector<SSENL<16>*> psse_tbl;
    LogMixerBias *pmix_ctx=nullptr;
};

#endif

