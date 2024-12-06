#include "pred.h"
#include <cassert>

Predictor::Predictor(const tparam &p)
:p(p),nA(p.nA),nB(p.nB),nM0(p.nM0),nS0(p.nS0),nS1(p.nS1),
ols{OLS(nA+nM0,p.k,p.lambda0,p.ols_nu0,p.beta_sum0,p.beta_pow0,p.beta_add0),
OLS(nB+nS0+nS1,p.k,p.lambda1,p.ols_nu1,p.beta_sum1,p.beta_pow1,p.beta_add1)},
lms{LMSCascade(p.vn0,p.vmu0,p.vmudecay0,p.vpowdecay0,p.mu_mix0,p.mu_mix_beta0),
LMSCascade(p.vn1,p.vmu1,p.vmudecay1,p.vpowdecay1,p.mu_mix1,p.mu_mix_beta1)},
be{BiasEstimator(p.bias_mu0,p.bias_scale),
   BiasEstimator(p.bias_mu1,p.bias_scale)}
{
}

void Predictor::fillbuf_ch0(const int32_t *src0,int idx0,const int32_t *src1,int idx1)
{
  vec1D &buf=ols[0].x;
  int bp=0;
  for (int i=idx0-nA;i<idx0;i++) buf[bp++]=(i>=0)?src0[i]:0.0;
  for (int i=idx1-nM0;i<idx1;i++) buf[bp++]=(i>=0)?src1[i]:0.0;
}

void Predictor::fillbuf_ch1(const int32_t *src0,const int32_t *src1,int idx1,int numsamples)
{
  vec1D &buf=ols[1].x;
  int bp=0;
  for (int i=idx1-nB;i<idx1;i++) buf[bp++]=(i>=0)?src1[i]:0.0;;
  for (int i=idx1-nS0;i<idx1+nS1;i++) buf[bp++]=(i>=0 && i<numsamples)?src0[i]:0.0;
}


double Predictor::predict(int ch)
{
  p_lpc[ch]=ols[ch].Predict();
  p_lms[ch]=lms[ch].Predict();
  return be[ch].Predict(p_lpc[ch]+p_lms[ch]);
}

void Predictor::update(int ch,double val)
{
  ols[ch].Update(val);
  lms[ch].Update(val-p_lpc[ch]);
  be[ch].Update(val);
}

