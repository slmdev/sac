#ifndef GLN2_H
#define GLN2_H

//Gated Linear Network
//mimics paq8 2-layer network with 1 hidden layer
//https://arxiv.org/abs/1910.01526

class GLN2 {
  public:
    //nctx = vector of hidden layer context sizes
    //nctx.size() = number of neurons
    GLN2(std::size_t nin,const std::vector<std::size_t>&nctx,float mconf=1.0f)
    :x(nin),
     nhidden(nctx.size()),
     fmix(nhidden),
     pred(nhidden),
     phidden(nhidden)
    {
      if (nhidden==0)
        throw std::invalid_argument("GLN2: hidden-layer size is zero");
      //nin mixer_weights at hidden[i][ctx]
      hidden.reserve(nhidden);
      for (auto ctx_size:nctx)
        hidden.emplace_back(ctx_size,LogitMixer(nin,mconf));
    }
    void SetConfidence(int i,float s)
    {
      for (auto &neuron:hidden)
        for (auto &mixer:neuron)
          mixer.SetConfidence(i,s);
    }
    int Predict(std::span<const int>inp,std::span<const int>hctx)
    {
      for (std::size_t i=0;i<inp.size();i++)
        x[i]=myDomain.Fwd(inp[i]);
      for (std::size_t i=0;i<nhidden;i++) {
        //select neuron based on context
        phidden[i]=&hidden[i][hctx[i]];
        pred[i] = phidden[i]->Predict(x);
      }
      return fmix.Predict(pred);
    }
    void Update(int bit,int hidden_rate,int final_rate)
    {
      for (auto *mixer:phidden)
        mixer->Update(x, bit,hidden_rate);
      fmix.Update(bit,final_rate);
    }
  protected:
    std::vector<int16_t>x;
    std::size_t nhidden;
    std::vector<std::vector<LogitMixer>> hidden;
    LogMixer fmix;
    std::vector<int>pred;
    std::vector<LogitMixer*> phidden;
};

#endif
