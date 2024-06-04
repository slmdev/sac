#ifndef PRED_H
#define PRED_H

#include "../pred/lms.h"
#include "../pred/lms_cascade.h"
#include "../pred/lpc.h"
#include "../pred/bias.h"

class Predictor {
  public:
    struct tparam {
      int nA,nB,nS0,nS1,k;
      std::vector <int>vn0,vn1;
      std::vector <double>vmu0,vmu1;
      std::vector <double>vmudecay0,vmudecay1;
      std::vector <double>vpowdecay0,vpowdecay1;
      double lambda0,lambda1,ols_nu0,ols_nu1,mu_mix0,mu_mix1,mu_mix_beta0,mu_mix_beta1;
      double beta_sum0,beta_pow0,beta_add0;
      double beta_sum1,beta_pow1,beta_add1;
      double mix_nu0,mix_nu1;
      double bias_mu0,bias_mu1;
      int bias_scale;
    };
    Predictor(const tparam &p);
    double Predict_stage0_ch0();
    void Update_stage0_ch0(double val);

    double Predict_stage0_ch1(const int32_t *ch_master,int nsample,int numsamples);
    void Update_stage0_ch1(double val);

    double PredictMaster();
    double PredictSlave(const int32_t *ch_master,int nsample,int numsamples);
    void UpdateMaster(double val);
    void UpdateSlave(double val);

    void FillSlaveHist(const int32_t *ch_master,int nsample,int numsamples,vec1D &buf);
    tparam p;
    int nA,nB,nS0,nS1;
    OLS<double>ols0,ols1;
    LMSCascade lms0,lms1;
    BiasEstimator be[2];
    vec1D hist0,hist1; //,tbuf;
    double p_lpc0,p_lpc1,p_lms0,p_lms1;
};

#endif // PRED_H
