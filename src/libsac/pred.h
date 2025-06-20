#ifndef PRED_H
#define PRED_H

#include "../pred/lms.h"
#include "../pred/lms_cascade.h"
#include "../pred/lpc.h"
#include "../pred/bias.h"

class Predictor {
  public:
    struct tparam {
      int nA,nB,nM0,nS0,nS1,k;
      std::vector <int>vn0,vn1;
      std::vector <double>vmu0,vmu1;
      std::vector <double>vmudecay0,vmudecay1;
      std::vector <double>vpowdecay0,vpowdecay1;
      double lambda0,lambda1,ols_nu0,ols_nu1,mu_mix0,mu_mix1,mu_mix_beta0,mu_mix_beta1;
      double beta_sum0,beta_pow0,beta_add0;
      double beta_sum1,beta_pow1,beta_add1;
      int ch_ref;
      double bias_mu0,bias_mu1;
      int bias_scale0,bias_scale1;
      int lm_n;
      double lm_alpha;
    };
    explicit Predictor(const tparam &p);

    double predict(int ch);
    void update(int ch,double val);

    void fillbuf_ch0(const int32_t *src0,int idx0,const int32_t *src1,int idx1);
    void fillbuf_ch1(const int32_t *src0,const int32_t *src1,int idx1,int numsamples);

    tparam p;
    int nA,nB,nM0,nS0,nS1;

    OLS ols[2];
    Cascade lms[2];
    BiasEstimator be[2];
    double p_lpc[2],p_lms[2];
};

#endif // PRED_H
