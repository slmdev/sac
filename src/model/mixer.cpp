#include "mixer.h"
#include "model.h"

using namespace BM;

LogMixer::LogMixer(int n)
:x(n),w(n),pd(0),n(n)
{
};

int LogMixer::Predict(const std::vector <int>&p)
{
  int64_t sum=0;
  for (int i=0;i<n;i++)
  {
    x[i]=myDomain.Fwd(p[i]);
    sum+=int64_t(w[i])*x[i];
  }
  sum=idiv_signed(sum,WBITS);
  pd=std::clamp(myDomain.Inv(sum),1,PSCALEm);
  return pd;
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

