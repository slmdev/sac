#include "pred.h"

Predictor::Predictor(const tparam &p)
:p(p),nA(p.nA),nS0(p.nS0),nS1(p.nS1),
ols0(nA,p.k,p.lambda0,p.ols_nu0),
ols1(nA+nS0+nS1+1,p.k,p.lambda1,p.ols_nu1),
lms0(p.vn0,p.vmu0,p.vmudecay0,p.vpowdecay0,p.mu_mix0,p.mu_mix_beta0),
lms1(p.vn1,p.vmu1,p.vmudecay1,p.vpowdecay1,p.mu_mix1,p.mu_mix_beta1),
be0(0.001,0.95),be1(0.001,0.95),
hist0(nA),hist1(nA),tbuf(nA+nS0+nS1+1)
{
  p_lpc0=p_lpc1=p_lms0=p_lms1=0.;
}

double Predictor::PredictMaster()
{
  p_lpc0=ols0.Predict(hist0);
  p_lms0=lms0.Predict();
  return be0.Predict(p_lpc0+p_lms0);
}

double Predictor::PredictSlave(const int32_t *ch_master,int nsample,int numsamples)
{
  FillSlaveHist(ch_master,nsample,numsamples,tbuf);
  p_lpc1=ols1.Predict(tbuf);
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
  //mix1.Update(val);
}

void Predictor::FillSlaveHist(const int32_t *ch_master,int nsample,int numsamples,vec1D &buf)
{
  std::fill(begin(buf),end(buf),0.);
  for (int k=0;k<nA;k++) buf[k]=hist1[k];
  int bp=nA;
  for (int k=-nS0;k<=nS1;k++) {
    int idx=nsample+k;
    if (idx>=0 && idx<numsamples) {
      buf[bp]=ch_master[idx];
    }
    bp++;
  }
}

