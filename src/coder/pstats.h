#ifndef PSTATS_H
#define PSTATS_H

#include <cstdint>
#include <vector>
#include <cassert>

class ExpMean {
public:
  explicit ExpMean(double alpha,int numsamples)
  :alpha_(alpha),count_(numsamples),
   right_sum_(numsamples),inv_weight_(numsamples)
  {
    assert(numsamples>0 && alpha>0 && alpha<1.0);
    //calc inv_weight inplace
    //pass1: inv_weight holds right_weight[i]
    //right_weight[i]=alpha+alpha^2 + ... + alpha^(count_-1-i)
    inv_weight_[count_-1]=0.0;
    for (int i=count_-2;i>=0;--i) {
      inv_weight_[i]=alpha_*(1.0 + inv_weight_[i+1]);
    }
    //pass2: fold in left_weight[i]=alpha + ... + alpha^i
    double left_weight_=0;
    for (int i=0;i<count_;++i) {
      inv_weight_[i]=1.0/(left_weight_ + 1.0 + inv_weight_[i]);
      left_weight_ = alpha_ * (left_weight_ + 1.0);
    }
  }

  void Begin(const int32_t* buf,int bpn)
  {
    assert(bpn >= 0 && bpn < 32);
    buf_ = buf;
    pos_ = 0;
    left_sum_ = 0.0;

    const uint32_t plane = uint32_t{1} << bpn;
    known_mask_ = ~(plane - 1);
    high_mask_ = ~(plane | (plane - 1));

    right_sum_[count_-1]=0.0;
    for (int i = count_-2;i >= 0; --i)
      right_sum_[i] = alpha_*(High(i+1)+right_sum_[i+1]);
  }

  double Mean() const
  {
    assert(pos_>=0 && pos_<count_);
    return (left_sum_+High(pos_)+right_sum_[pos_])*inv_weight_[pos_];
  }

  void Advance()
  {
    assert(pos_<(count_-1));
    left_sum_ = alpha_ * (left_sum_ + Known(pos_));
    ++pos_;
  }

private:
  inline double Known(int i) const {
    return static_cast<uint32_t>(buf_[i]) & known_mask_;
  }

  inline double High(int i) const {
    return static_cast<uint32_t>(buf_[i]) & high_mask_;
  }

  double alpha_;
  int count_;
  std::vector<double> right_sum_, inv_weight_;
  const int32_t* buf_ = nullptr;
  int pos_ = 0;
  uint32_t known_mask_ = 0, high_mask_ = 0;
  double left_sum_ = 0;
};

class SigStats {
  public:
    struct St {
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int n4 = 0;
    };
    SigStats(int radius,int max_age=2)
    :radius_(radius),max_age_(max_age)
    {

    }
    //n1=number of significant samples
    //n2=number of significant samples with age>0
    //n3=sum of capped_age
    //n4=sum of distance*capped_age
    St Counts() const
    {
      return {
        left_.n1 + right_.n1,
        left_.n2 + right_.n2,
        left_.n3 + right_.n3,
        left_.n4 + right_.n4
      };
    }
    //max of activity sum n4
    int get_max_acc() {
      return max_age_*radius_*(radius_+1);
    }
    void Begin(const std::vector<int>&msb,int bpn)
    {
      msb_= msb.data();
      pos_= 0;
      count_=msb.size();
      bpn_= bpn;
      left_={};
      right_={};
      const int k=std::min(radius_,count_-1);
      for (int j=1;j<=k;++j)
        Change(right_,j,+1);
    }
    void Advance()
    {
      assert(pos_ + 1 < count_);
      // Existing left samples become one position farther away.
      left_.n4 += left_.n3;
      // Existing right samples become one position closer.
      right_.n4 -= right_.n3;

      ++pos_;

      // Remove the sample now beyond the left boundary.
      Change(left_,pos_-radius_-1, -1);
      // The processed center enters the left window at distance 1.
      Change(left_,pos_-1 ,+1);
      // The next center leaves the right window at distance zero.
      Change(right_,pos_,-1);
      // A new sample enters at the far right.
      Change(right_,pos_+radius_,+1);
    }
  protected:
    void Change(St& st, int index, int sign)
    {
      if (index < 0||index >= count_)
          return;
      const int x = msb_[index];
      if (x == 0) return;
      assert(x >= bpn_);
      st.n1 += sign;
      const int d = std::min(x-bpn_,max_age_);
      if (d==0) return;
      const int distance=std::abs(index-pos_);
      st.n2 += sign;
      st.n3 += sign*d;
      st.n4 += sign*d*distance;
    }
    int radius_,max_age_;
    const int* msb_ = nullptr;
    int count_=0;
    int bpn_=0;
    int pos_=0;
    St left_,right_;
};

#endif

