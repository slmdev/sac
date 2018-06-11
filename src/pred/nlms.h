#ifndef NLMS_H
#define NLMS_H

#include "../common/histbuf.h"
#include "blend.h"
#include "lpc.h"

#define APRIORI_UPDATE  //updates energy before updating weighting coefs
//#define CLAMP_PRED // clamp the prediction to values already seen
#define NLMS_EPS (0.0001)
#define MOMENTUM_ALPHA (0.0) //adds a momentum to the velocity (typically 0.5-0.9)

class BlendLM {
  public:
      BlendLM(double lambda,int n,int kmax)
      :myLM(lambda,n,kmax),n(n)
      {

      }
      double Predict(const std::vector<double> &p)
      {
        for (int i=0;i<n;i++) myLM.x[i]=p[i];
        return myLM.Predict();
      }
      void Update(double val)
      {
        myLM.Update(val);
      }
  private:
    LM myLM;
    int n;
};

inline float rsqrt(float __x)
{
    float reciprocal;
    __asm__ __volatile__ (
        "movss %1, %%xmm0\n"
        "rsqrtss %%xmm0, %%xmm1\n"
        "movss %%xmm1, %0\n"
        :"=m"(reciprocal)
        :"m"(__x)
        :"xmm0", "xmm1"
    );

    //
    // One iteration of Newton's method. Doubles the running time, but adds
    // some MUCH needed accuracy to the inline SSE asm instruction.
    //
    //return (0.5f * (reciprocal + 1.0f/(__x * reciprocal)));
    //
    // If we need this to run 2-3 times as fast, then use the line below.
    // This causes the function to lose a significant amount of accuracy.
    //
    return reciprocal;
    //
}

// ADAGRAD method with adaptive subgradient
// the inverse squareroot makes it horrible slow
class LMSADA {
  public:
      LMSADA(double mu,int order)
      :w(order),eg(order),hist(order),order(order),pred(0),mu(mu)
      {
      }
      int32_t Predict() // calculate prediction
      {
        double sum=0.;
        for (int i=0;i<order;i++) sum+=w[i]*hist[i];
        pred=static_cast<int32_t>(floor(sum+0.5));
        return pred;
      }
      void Update(int32_t val)
      {
        const double err=val-pred; // prediction error
        const double rho=0.95;
        for (int i=0;i<order;i++) {
            double const grad=err*hist[i]; // gradient
            eg[i]=rho*eg[i]+(1.0-rho)*grad*grad; //accumulate gradients
            //w[i]+=(mu*grad)/sqrt(eg[i]+1.0);// update weights
            w[i]+=(mu*grad*rsqrt(eg[i]+1.0));// update weights
            //w[i]+=(mu*grad/(hist[i]*hist[i]+0.001));// update weights
            //w[i]+=(mu*grad/(eg[i]+1.));// update weights
        }
        hist.PushBack(val);
      }
      ~LMSADA(){};
  private:
    std::vector<double> w,eg;
    HistBuffer <double>hist;
    int32_t order,pred;
    double mu;
};

class LMSADA2 {
inline double sgn(double v){
      if (v<0) return -1;
      if (v>0) return 1;
      return 0;
    }
  public:
      LMSADA2(int S,int D,double mu1,double mu2)
      :w(S+D),eg(S+D),hist(S+D),pred(0),mu1(mu1),mu2(mu2),S(S),D(D)
      {
      }
      double Predict(double val1) // calculate prediction
      {
        for (int i=(S+D)-1;i>S;i--) hist[i]=hist[i-1];hist[S]=val1; // update history for right channel
        pred=0.;
        for (int i=0;i<S+D;i++) pred+=w[i]*hist[i];
        return pred;
      }
      void Update(double val)
      {
        const double err=val-pred; // prediction error
        const double rho=0.95;
        double mu=mu1;
        for (int i=0;i<S;i++) {
            double const grad=err*hist[i]; // gradient
            eg[i]=rho*eg[i]+(1.0-rho)*(grad*grad);
            w[i]+=(mu*grad*rsqrt(eg[i]+0.001));// update weights
        }
        mu=mu2;
        for (int i=S;i<S+D;i++) {
            double const grad=err*hist[i]; // gradient
            eg[i]=rho*eg[i]+(1.0-rho)*(grad*grad);
            w[i]+=(mu*grad*rsqrt(eg[i]+0.001));// update weights
        }
        for (int i=S-1;i>0;i--) hist[i]=hist[i-1];hist[0]=val; //update history for left channel
      }
      ~LMSADA2(){};
  private:
    std::vector<double> w,eg;
    std::vector<double> hist;
    double pred;
    double mu1,mu2;
    int S,D;
};

// classical NLMS algorithm
class NLMS {
  public:
      NLMS(double mu,int order)
      :w(order),hist(order),order(order),pred(0),spow(0.0),mu(mu/(double)order)
      {
      }
      int32_t Predict()
      {
        pred=0.; // caluclate prediction
        for (int i=0;i<order;i++) pred+=w[i]*hist[i];
        return pred;
      }
      void Update(double val) {

        double err=val-pred; // calculate prediction error

        spow-=(hist[order-1]*hist[order-1]); // update energy in the tapping window
        spow+=(val*val);
        const double  wf=mu*err/(0.001+spow); //constant weight factor
        for (int i=0;i<order;i++) w[i]+=wf*hist[i];
        hist.PushBack(val);
      }
      ~NLMS(){};
  private:
    std::vector<double> w;
    HistBuffer <double>hist;
    int32_t order;
    double pred,spow,mu;
};


class NLMS2 {
    inline double sgn(double v){
      if (v<0) return -1;
      if (v>0) return 1;
      return 0;
    }
  public:
      NLMS2(double mu,int S,int D):w(S+D),dw(S+D),hist(S+D),mu(mu),S(S),D(D) {
        pred=0;
        spow=0.0;
      }
      int32_t Predict(int32_t val1)
      {
        UpdateEnergy( (S+D)-1,val1);
        for (int i=(S+D)-1;i>S;i--) hist[i]=hist[i-1];hist[S]=val1; // update history for right channel

        double sum=0.; // prediction d(n)=w(n-1)^T*x(n)
        for (int i=0;i<S+D;i++) sum+=w[i]*hist[i];
        pred=floor(sum+0.5);
        return pred;
      }
      void Update(int32_t val) {
        #ifdef APRIORI_UPDATE
          UpdateEnergy(S-1,val);
        #endif

        const int32_t err=val-pred; // prediction error

        spow=0.0;
        double beta=1.;
        for (int i=0;i<S+D;i++) {
         spow+=beta*(hist[i]*hist[i]);
         beta*=0.98;
        }
        // update weight vector w(n+1)=w(n)+(mu*err*x(n)))/(x(n)^T*x(n))
        const double  wf=mu*sgn(err)/(spow+NLMS_EPS); //constant weight factor
        for (int i=0;i<S+D;i++) {
                dw[i]=MOMENTUM_ALPHA*dw[i]+wf*hist[i];
                w[i]+=dw[i];

        }

        #ifndef APRIORI_UPDATE
          UpdateEnergy(S-1,val);
        #endif
        for (int i=S-1;i>0;i--) hist[i]=hist[i-1];hist[0]=val; //update history for left channel
      }
      ~NLMS2(){};
  private:
    void UpdateEnergy(int pos,int32_t val) {
     spow-=(hist[pos]*hist[pos]); // update energy in the tapping window
     spow+=(val*val);
    }
    std::vector<double> w,dw,hist;
    double spow,mu;
    int S,D,pred;
};

#endif // NLMS_H
