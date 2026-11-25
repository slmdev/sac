#include "mixer.h"
#include "model.h"

int LogMixer::Predict(std::span<const int>p)
{
  for (int i=0;i<n;i++)
    x[i]=myDomain.Fwd(p[i]);

  const int logit=dot(x,w);
  return (pd=myDomain.Inv(logit));
}


int LogMixerBias::Predict(int p1,std::span<const int>c)
{
  int base = myDomain.Fwd(p1);
  for (int i=0;i<n;i++) {
    int z=myDomain.Fwd(c[i]);
    x[i]=(z-base);
  }
  const int logit=dot(x,w);
  pd=myDomain.Inv(base+logit);
  return pd;
}
