#ifndef RBF_H
#define RBF_H

#include "lpc.h"

class Map01 {
  public:
      Map01(double vmin,double vmax)
      :vmin(vmin),vmax(vmax),vrange(vmax-vmin)
      {

      }
      double Map(double val)
      {
         return (val-vmin)/vrange;
      }
      double Unmap(double val) {
        return val*vrange+vmin;
      }
  private:
    double vmin,vmax,vrange;
};

// classical k-means
class KMeans {
  public:
      KMeans(int k,int n)
      :h(k,std::vector<double>(n)),hm(k,std::vector<double>(n)),nh(k),sigma(k),n(n),k(k)
      {
        srand(123456789L);
      }
      int rand_int(int imin,int imax) {
        return imin + ( std::rand() % ( imax - imin + 1 ) );
      }
      void PrintVec(const std::vector <double>&vec)
      {
         for(auto i : vec) std::cout << i << " ";
         std::cout << std::endl;
      }
      void InitCenters(const std::vector <std::vector<double>>&data)
      {
         for (int i=0;i<k;i++) {
            int r=rand_int(0,data.size()-1);
            h[i]=data[r];
         }
      }
      void RunKMeans(const std::vector <std::vector <double>>&data)
      {
        for (int i=0;i<k;i++) {nh[i]=1;sigma[i]=0.;};

        for (size_t i=0;i<data.size();i++) {
           const std::vector <double>&val=data[i];

           double mind=MathUtils::L2Dist(val,h[0]); // find cluster with minimum distance
           int minc=0;
           for (int c=1;c<k;c++) {
             double d=MathUtils::L2Dist(val,h[c]);
             if (d<mind) {mind=d;minc=c;};
           }
           nh[minc]++;
           double alpha=1.0/(double)nh[minc]; // update cluster
           //double alpha=0.001;
           std::vector <double>&bmu=h[minc];
           for (int c=0;c<n;c++) bmu[c]=(1.0-alpha)*bmu[c]+alpha*val[c];
           sigma[minc]=(1.0-alpha)*sigma[minc]+alpha*mind;
        }
      }
      double CalcMeanCenterDist()
      {
        double sum=0.;
        for (int i=0;i<k;i++) sum+=MathUtils::L2Dist(h[i],hm[i]);
        sum/=(double)k;
        return sum;
      }
      int CalcCenters(const std::vector <std::vector<double>>&data)
      {
        InitCenters(data);
        hm=h;
        int num_iterations=0;
        while (1) {
          RunKMeans(data);
          num_iterations++;
          if (CalcMeanCenterDist()<EPS || num_iterations>100) break;
          hm=h;
        }
        return num_iterations;
        //cout << "iterations: " << num_iterations << endl;
        //for (auto &v:h) PrintVec(v);
        //cout << endl;
      }
    std::vector<std::vector<double>> h;
  private:
    std::vector<std::vector<double>> hm;
    std::vector <int>nh;
    std::vector <double>sigma;
    int n,k;
};

class BatchRBF {
  public:
    BatchRBF(int num_centers,int num_memory,int batch_size)
    :myLM(0.999,num_centers+2),
     myKMeans(num_centers,num_memory),
     myMap(-(1<<15),1<<15),
     data(batch_size,std::vector<double>(num_memory)),
     hist(num_memory),buf(batch_size),sigma(num_centers,1.0),
     ncenters(num_centers),nmemory(num_memory),bsize(batch_size)
    {
      k=0;centers_avail=false;
    }
    void FillLM(const std::vector<double> &v)
    {
      for (int i=0;i<ncenters;i++)
        myLM.x[i]=Radial(MathUtils::L2Dist(v,myKMeans.h[i]),sigma[i]);
      myLM.x[ncenters]=v[0];
      myLM.x[ncenters+1]=v[1];
      //myLM.x[ncenters+2]=1;
    }
    double Predict(){
      FillLM(hist);
      double p=myLM.Predict();
      return p;
    }
    double Predict(const std::vector<double> &v) {
       FillLM(v);
       return myLM.Predict();
    }
    double Radial(double d,double s)
    {
       return exp(-(d*d)/(2*s*s));
       //return d*d*log(d);
    }
    void CalcSigma() {
      double avgd=0.;
      for (int j=0;j<ncenters;j++) {
        double mind1=std::numeric_limits<double>::max();
        double mind2=std::numeric_limits<double>::max();
        for (int i=0;i<ncenters;i++) {
          if (i!=j) {
            double d=MathUtils::L2Dist(myKMeans.h[i],myKMeans.h[j]);
            if (d<mind1) mind1=d;
            else if (d<mind2) mind2=d;
            avgd+=d;
          }
        }
        sigma[j]=(mind1+mind2)/2.;
      }
      avgd/=double(ncenters*(ncenters-1));
      for (int i=0;i<ncenters;i++) sigma[i]=avgd;
    }
    void CalcRBFWeights() {
      //myLM.Init();
      for (int i=0;i<bsize;i++) {
        FillLM(data[i]);
        myLM.Update(buf[i]);
      }
      myLM.Solve();
    }
    void TestRBF() {
      double ivar=0.;
      double evar=0.;
      for (int i=0;i<bsize;i++) {
        double p=Predict(data[i]);
        double e=buf[i]-p;
        ivar+=buf[i]*buf[i];
        evar+=e*e;
      }
      if (bsize) {
        ivar/=(double)bsize;
        evar/=(double)bsize;
        double gain=10.0*log(ivar/evar)/log(10.);
        std::cout << "gain: " << gain << " dB\n";
      }
    }
    void Update(double val)
    {
       data[k]=hist;
       buf[k]=val;
       k++;
       if (k==bsize) {
         k=0;
         myKMeans.CalcCenters(data);
         CalcSigma();
         //cout << "rbf-weights...\n";
         CalcRBFWeights();
         if (!centers_avail) centers_avail=true;
         //cout << "test rbf...\n";
         //TestRBF();
       }
       for (int i=nmemory-1;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
    }
  private:
    LMM myLM;
    KMeans myKMeans;
    Map01 myMap;
    std::vector <std::vector<double>>data;
    std::vector <double>hist,buf,sigma;
    int ncenters,nmemory,bsize,k;
    bool centers_avail;
};

class RBFFuzzy {
  public:
      RBFFuzzy(int num_centers,int num_memory)
      :myLM(0.998,num_centers*num_memory,1000),
       myMap(-(1<<15),1<<15),
       h(num_memory,std::vector<double>(num_centers,0)),
       hist(num_memory),histm(10),sigma(num_memory,1),
       ncenters(num_centers),nmemory(num_memory)
      {

      }
      void UpdateCenters(double nu)
      {
         for (int i=0;i<nmemory;i++)  {
            double x=hist[i];
            double mind=std::numeric_limits<double>::max();
            double minc=-1;
            for (int j=0;j<ncenters;j++){
                double d=fabs(x-h[i][j]);
                if (d<mind) {
                    mind=d;
                    minc=j;
                }
            }
            h[i][minc]=nu*h[i][minc]+(1.0-nu)*x;

            double avgd=0.;
            for (int k=0;k<ncenters;k++)
              for (int l=0;l<ncenters;l++) {
                if (l!=k) avgd+=fabs(h[i][k]-h[i][l]);
              }
            avgd/=double(ncenters*(ncenters-1));
            sigma[i]=avgd;
         }
      }
      double Radial(double d,double s)
      {
        return exp(-d*d/(2*s*s));
      }
      double Predict() {
        int k=0;
        for (int i=0;i<nmemory;i++) {
           for (int j=0;j<ncenters;j++) {
              myLM.x[k++]=(histm[i]+h[i][j])*Radial(fabs(hist[i]-h[i][j]),10000);
           }
        }
        return myLM.Predict();
      }
      void Update(double val)
      {
       UpdateCenters(0.998);
       myLM.Update(val);
       for (int i=9;i>0;i--) histm[i]=histm[i-1];histm[0]=val;
       hist[0]=histm[0]-histm[1];
       //for (int i=nmemory-1;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
      }
  private:
    LM myLM;
    Map01 myMap;
    std::vector <std::vector<double>> h;
    std::vector <double>hist,histm,sigma;
    int ncenters,nmemory;
};

#endif // RBF_H
