#include "mixer.h"
#include "model.h"

using namespace BM;

LogMixer::LogMixer(int n,InitType type)
:x(n),w(n),pd(0),n(n)
{
  if (n && type==InitType::Uniform)
    for (auto &weight:w)
      weight=WSCALE/n;
};

int LogMixer::Predict(const std::vector<int>&p)
{
  for (int i=0;i<n;i++)
    x[i]=myDomain.Fwd(p[i]);

  const int pl=dot(x,w);
  return (pd=myDomain.Inv(pl));
}

void LogMixer::Update(int bit,int rate)
{
  int err=(bit<<PBITS)-pd;
  for (int i=0;i<n;i++)
  {
    int dw=idiv_signed(int64_t(x[i])*err*rate,USHIFT);
    w[i]=std::clamp(w[i]+dw,-WRANGE,WRANGE-1);
  }
};

