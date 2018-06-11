#ifndef LPC_H
#define LPC_H

#include "../math/cholesky.h"
#include "../math/cov.h"
#include "blend.h"


class LPC {
  public:
      LPC(double alpha,int order,int k)
      :cov(alpha,order),chol(order),sol(order),order(order),k(k),nk(0)
      {

      }
      int32_t Predict() {
        double sum=0.;
        for (int i=0;i<order;i++) sum+=sol[i]*cov.hist[i];
        return floor(sum+0.5);
      }
      void Update(int32_t val) {
        cov.UpdateCov();
        cov.UpdateB(val);
        cov.UpdateHist(val,0,order-1);
        nk++;
        if (nk>=k) {
          nk=0;
          if (chol.Factor(cov.c)) {
            std::cout << "matrix indefinit" << std::endl;
          } else chol.Solve(cov.b,sol);
        }
      }
  //private:
    AutoCov cov;
    Cholesky chol;
    std::vector<double> sol;
    int order,k,nk;
};

// general linear model of the form p=w1*p1+w2*p2+...+wn*pn
class LM {
  public:
    LM(double lambda,int n,int kmax)
    :x(n),w(n),b(n),cov(n),chol(n),lambda(lambda),n(n),kmax(kmax)
    {
      k=0;
      for (int i=0;i<n;i++) cov[i][i]=1.0;
    }
    double Predict() {
      double sum=0.;
      for (int i=0;i<n;i++) sum+=w[i]*x[i];
      return sum;
    }
    void Update(double val)
    {
      int i,j;
      for (j=0;j<n;j++)
        for (i=0;i<n;i++) cov[j][i]=lambda*cov[j][i]+(x[j]*x[i]);

      for (i=0;i<n;i++) b[i]=lambda*b[i]+(x[i]*val);

      k++;
      if (k>=kmax) {
        if (chol.Factor(cov)) {
          std::cout << "matrix indefinit" << std::endl;
        } else chol.Solve(b,w);
        k=0;
      }
    }
    std::vector<double> x;
  private:
    std::vector<double> w,b;
    Matrix cov;
    Cholesky chol;
    double lambda;
    int n,kmax,k;
};

// general linear model of the form p=w1*p1+w2*p2+...+wn*pn
class LMM {
  public:
    LMM(double lambda,int n)
    :x(n),w(n),b(n),cov(n),chol(n),lambda(lambda),n(n)
    {
      Init();
    }
    void Init() {
      for (int j=0;j<n;j++)
      for (int i=0;i<n;i++) {
        if (i==j) cov[j][i]=1.0;
        else cov[j][i]=0.;
      }
      for (int i=0;i<n;i++) {w[i]=0;b[i]=0;};
    }
    double Predict() {
      double sum=0.;
      for (int i=0;i<n;i++) sum+=w[i]*x[i];
      return sum;
    }
    void Update(double val)
    {
      int i,j;
      for (j=0;j<n;j++)
        for (i=0;i<n;i++) cov[j][i]=lambda*cov[j][i]+(x[j]*x[i]);

      for (i=0;i<n;i++) b[i]=lambda*b[i]+(x[i]*val);
    }
    void Solve() {
        if (chol.Factor(cov)) {
          std::cout << "matrix indefinit" << std::endl;
          for (int i=0;i<n;i++) w[i]=0.;
        } else chol.Solve(b,w);
    }
    std::vector<double> x,w;
  private:
    std::vector<double> b;
    Matrix cov;
    Cholesky chol;
    double lambda;
    int n;
};

class LM1 {
  public:
      LM1(double lambda,int n,int kmax)
      :myLM(lambda,n,kmax),hist(n),n(n)
      {

      }
      int32_t Predict() {
        for (int i=0;i<n;i++) myLM.x[i]=hist[i];
        double pm=myLM.Predict();
        return round(pm);
      }
      void Update(int32_t val)
      {
        myLM.Update(val);
        for (int i=n;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
      }
  private:
    LM myLM;
    std::vector <int32_t>hist;
    double mean;
    int n;
};

// Volterra Filter of order-2 with history n
class LMVolt2 {
  public:
      LMVolt2(double lambda,int n,int kmax)
      :n(n),myLM(lambda,(n*(n+1))/2,kmax),hist(n)
      {
      }
      double Predict() {
        int k=0;
        for (int j=0;j<n;j++)
          for (int i=0;i<=j;i++) myLM.x[k++]=(hist[i]*hist[j]);
        return myLM.Predict();
      }
      void Update(double val) {
        myLM.Update(val);
        for (int i=n;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
      }
  protected:
    int n;
    LM myLM;
    std::vector<double> hist;
};

//#define REORDER_HIST

class LPC3 {
  public:
int getIdxA(int i) {
  return 3*i/2;
}

int getIdxB(int i) {
  return 3*i+2;
}
      LPC3(double alpha,int S,int D,int k,double regp=0.)
      :cov(alpha,S+D),chol(S+D),sol(S+D),ssum(S+D),blend(S+D,0.98),S(S),D(D),k(k),nk(0)
      {
        chol.regp=regp;
      }
      double Predict(double val1) {
        #ifdef REORDER_HIST
        if (D) {
          for (int i=D-1;i>0;i--) cov.hist[getIdxB(i)]=cov.hist[getIdxB(i-1)];
          cov.hist[getIdxB(0)]=val1;
        }
        #else
         if (D) cov.UpdateHist(val1,S,(S+D)-1);
        #endif
        for (int order=0;order<S+D;order++) {
          double sum=0.;
          for (int i=0;i<=order;i++) sum+=chol.sol[order][i]*cov.hist[i];
          ssum[order]=sum;
        }
        return blend.Predict(ssum);
      }
      void Update(double val0) {
        blend.Update(val0);
        cov.UpdateCov();
        cov.UpdateB(val0);
        #ifdef REORDER_HIST
        for (int i=S-1;i>0;i--) cov.hist[getIdxA(i)]=cov.hist[getIdxA(i-1)];
        cov.hist[getIdxA(0)]=val0;
        #else
        cov.UpdateHist(val0,0,S-1);
        #endif

        nk++;
        if (nk>=k) {
          nk=0;
          if (chol.Factor(cov.c)) {
            std::cout << "matrix indefinit" << std::endl;
          } else chol.Solve(cov.b);
        }
      }
      //~LPC2(){cout << "lpc2 destruct\n";};
  private:
    AutoCov cov;
    Cholesky chol;
    std::vector<double> sol,ssum;
    BlendErr blend;
    int S,D,k,nk;
};

class VF1 {
  public:
    VF1(double alpha,double beta)
    :ecor(alpha),epow(alpha),p(0.998),lerr(0.),beta(beta)
    {
      lmax=0.9999;lmin=0.99;
    }
    double Get() {
      return p;
    }
    void Update(double err) {
      double r=ecor.Get()/(epow.Get()+EPS);
      double t=exp(-beta*r*r);
      p=lmin+(lmax-lmin)*t;

      ecor.Update(err*lerr);
      epow.Update(err*err);
      lerr=err;
    }
  private:
    RunExp ecor,epow;
    double lmax,lmin,p,lerr,beta;
};

class LPC2 {
  public:
      LPC2(double alpha,int S,int D,int k,double regp=0.)
      :cov(alpha,S+D),chol(S+D),sol(S+D),S(S),D(D),k(k),nk(0)
      {
         chol.regp=regp;
      }
      double Predict(double val1) {
        if (D) cov.UpdateHist(val1,S,(S+D)-1);

        pred=0.;
        for (int i=0;i<S+D;i++) pred+=sol[i]*cov.hist[i];
        return pred;
      }
      void Update(double val) {
        cov.UpdateCov();
        cov.UpdateB(val);
        cov.UpdateHist(val,0,S-1);

        nk++;
        if (nk>=k) {
          nk=0;
          if (chol.Factor(cov.c)) {
            std::cout << "matrix indefinit" << std::endl;
          } else chol.Solve(cov.b,sol);
        }
      }
      //~LPC2(){cout << "lpc2 destruct\n";};
  private:
    AutoCov cov;
    Cholesky chol;
    std::vector<double> sol;
    int S,D,k,nk;
    double pred;
};

#endif // LPC_H
