#include "pred.h"

Predictor::Predictor(const tparam &p)
:p(p),nA(p.nA),nB(p.nB),nS0(p.nS0),nS1(p.nS1),
ols0(nA,p.k,p.lambda0,p.ols_nu0,p.beta_sum0,p.beta_pow0,p.beta_add0),
ols1(nB+nS0+nS1+1,p.k,p.lambda1,p.ols_nu1,p.beta_sum0,p.beta_pow0,p.beta_add0),
lms0(p.vn0,p.vmu0,p.vmudecay0,p.vpowdecay0,p.mu_mix0,p.mu_mix_beta0,p.mix_nu0),
lms1(p.vn1,p.vmu1,p.vmudecay1,p.vpowdecay1,p.mu_mix1,p.mu_mix_beta1,p.mix_nu1),
be0(p.bias_mu,p.bias_rescale),be1(p.bias_mu,p.bias_rescale),
hist0(nA),hist1(nB)
{
  p_lpc0=p_lpc1=p_lms0=p_lms1=0.;
}

double Predictor::PredictMaster()
{
  for (int i=0;i<nA;i++) ols0.x[i]=hist0[i];
  p_lpc0=ols0.Predict();
  p_lms0=lms0.Predict();
  return be0.Predict(p_lpc0+p_lms0);
}

double Predictor::PredictSlave(const int32_t *ch_master,int nsample,int numsamples)
{
  FillSlaveHist(ch_master,nsample,numsamples,ols1.x);
  p_lpc1=ols1.Predict();
  p_lms1=lms1.Predict();
  return be1.Predict(p_lpc1+p_lms1);
}

void Predictor::UpdateMaster(double val)
{
  ols0.Update(val);
  lms0.Update(val-p_lpc0);
  miscUtils::RollBack(hist0,val);
  be0.Update(val);
}

void Predictor::UpdateSlave(double val)
{
  ols1.Update(val);
  lms1.Update(val-p_lpc1);
  miscUtils::RollBack(hist1,val);
  be1.Update(val);
}

void Predictor::FillSlaveHist(const int32_t *ch_master,int nsample,int numsamples,vec1D &buf)
{
  std::fill(begin(buf),end(buf),0.);
  for (int k=0;k<nB;k++) buf[k]=hist1[k];
  int bp=nB;
  for (int k=-nS0;k<=nS1;k++) {
    int idx=nsample+k;
    if (idx>=0 && idx<numsamples) {
      buf[bp]=ch_master[idx];
    }
    bp++;
  }
}

