#include "profile.h"

int SacProfile::LoadBaseProfile()
{
  const int mo_lpc=32; // maximum ols order
  const int wbits_lms=13; // maximum nlms order 2^wbits_lms
  SacProfile &profile=*this;

  profile.Init(53);

  profile.Set(0,0.99,0.9999,0.998);
  profile.Set(1,0.001,10.0,0.1);

  profile.Set(2,0.001,1.0,0.1);//mu0
  profile.Set(3,0.001,1.0,0.12);//mu1
  profile.Set(4,0.001,1.0,0.06);//mu2
  profile.Set(5,0.001,1.0,0.04);//mu3

  profile.Set(6,0.98,1,1.0); // mu-decay
  profile.Set(7,0.0,1.0,0.8);   // pow-decay
  profile.Set(8,0.0,1.0,0.8);   // pow-decay

  profile.Set(10,0.0001,0.008,0.002);//mu-mix
  profile.Set(11,0.8,0.9999,0.95);//mu-mix-beta

  profile.Set(12,0.99,0.9999,0.998);
  profile.Set(13,0.001,10.0,0.1);

  profile.Set(14,0.001,1.0,0.1);//mu0
  profile.Set(15,0.001,1.0,0.12);//mu1
  profile.Set(16,0.001,1.0,0.06);//mu2
  profile.Set(17,0.001,1.0,0.04);//mu3

  profile.Set(18,0.98,1,1.0); // mu-decay
  profile.Set(19,0.0,1.0,0.8);   // pow-decay
  profile.Set(20,0.0,1.0,0.8);   // pow-decay
  profile.Set(21,0.0,1.0,0.8);   // pow-decay
  profile.Set(22,0.0001,0.008,0.002);//mu-mix
  profile.Set(23,0.8,0.9999,0.95);//mu-mix-beta*/

  profile.Set(24,4,mo_lpc,16);//nA
  profile.Set(25,4,mo_lpc,16);//nB
  profile.Set(26,0,mo_lpc,8);//nS0
  profile.Set(27,-mo_lpc,mo_lpc,8);//nS1
  profile.Set(9,0,mo_lpc,0); //nM0

  profile.Set(28,256,1<<wbits_lms,1280);
  profile.Set(29,32,1<<(wbits_lms-1),256);
  profile.Set(30,4,1<<(wbits_lms-2),32);

  profile.Set(31,256,1<<wbits_lms,1280);
  profile.Set(32,32,1<<(wbits_lms-1),256);
  profile.Set(33,4,1<<(wbits_lms-2),32);

  profile.Set(34,0,1,0.6);
  profile.Set(35,0.1,2,0.8);
  profile.Set(36,0,10,2);

  profile.Set(37,2,1<<(wbits_lms-3),4); //stage 4
  profile.Set(38,2,1<<(wbits_lms-3),4);

  profile.Set(39,0.98,1,1.0); // mu-decay
  profile.Set(40,0.98,1,1.0); // mu-decay

  profile.Set(41,1,10,4); //stage-5 lm
  profile.Set(42,0.1,10.0,5); // shape parameter gamma


  profile.Set(43,0.001,0.005,0.0015);//bc-mu0
  profile.Set(44,0.001,0.005,0.0015);//bc-mu1

  profile.Set(45,4,10,5);//bias scale in bits

  profile.Set(46,0.98,1,1.0); // mu_decay
  profile.Set(48,0.98,1,1.0); // mu_decay
  profile.Set(50,0.0,1.0,0.8); //pow_decay
  profile.Set(47,0.98,1,1.0); // mu_decay
  profile.Set(49,0.98,1,1.0); // mu_decay
  profile.Set(51,0.0,1.0,0.8); //pow_decay
  profile.Set(52,0.0,1.0,0.8); //pow_decay

  return profile.coefs.size();
}
