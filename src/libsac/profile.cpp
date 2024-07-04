#include "profile.h"

void LoadProfileNormal(SacProfile &profile)
{
    const int mo_lpc=32; // maximum ols order
    profile.Init(35,0);
    profile.Set(0,0.99,0.9999,0.998); // ols-lambda
    profile.Set(1,0.0001,10.0,0.001);//ols-nu
    profile.Set(2,0.00001,0.008,0.0008);//mu0
    profile.Set(3,0.00008,0.008,0.002);//mu1
    profile.Set(4,0.0001,0.05,0.008);//mu2
    profile.Set(5,0.0,1.0,0.8);   // pow-decay0
    profile.Set(6,0.0,1.0,0.8);   // pow-decay1
    profile.Set(7,0.0,1.0,0.8);   // pow-decay2
    profile.Set(8,0.98,1,1.0); // mu-decay0
    profile.Set(9,0.98,1,1.0); // mu-decay1
    profile.Set(10,0.98,1,1.0); // mu-decay2
    profile.Set(11,0.0001,0.008,0.002);//mu-mix
    profile.Set(12,0.8,0.9999,0.95); //mu-mix-beta

    profile.Set(13,0.99,0.9999,0.998); // ols-lambda
    profile.Set(14,0.001,10.0,0.001);//ols-nu
    profile.Set(15,0.00001,0.008,0.0008);//mu0
    profile.Set(16,0.00008,0.008,0.002);//mu1
    profile.Set(17,0.0001,0.05,0.008);//mu2
    profile.Set(18,0.0,1.0,0.8);   // pow-decay0
    profile.Set(19,0.0,1.0,0.8);   // pow-decay1
    profile.Set(20,0.0,1.0,0.8);   // pow-decay2
    profile.Set(21,0.98,1,1.0); // mu-decay0
    profile.Set(22,0.98,1,1.0); // mu-decay1
    profile.Set(23,0.98,1,1.0); // mu-decay2
    profile.Set(24,0.0001,0.008,0.002);//mu-mix
    profile.Set(25,0.8,0.9999,0.95); //mu-mix-beta

    profile.Set(26,4,mo_lpc,16);//ols-order
    profile.Set(27,4,mo_lpc,16);//ols-order
    profile.Set(28,0,mo_lpc,8);//ols-order
    profile.Set(29,0,mo_lpc,8);//ols-order

    profile.Set(30,0,1,0.6);//ols-order
    profile.Set(31,0.8,2,0.8);//ols-order
    profile.Set(32,0,10,2);//ols-order

    profile.Set(33,1,8,5);//bias_scale

    profile.Set(34,0.1,1.0,1.0);//mix_nu
}

void LoadProfileHigh(SacProfile &profile)
{
    const int mo_lpc=32; // maximum ols order

    profile.Init(53,1);

    profile.Set(0,0.99,0.9999,0.998);
    profile.Set(1,0.001,10.0,0.001);

    profile.Set(2,0.001,1.0,0.1);//mu0
    profile.Set(3,0.001,1.0,0.12);//mu1
    profile.Set(4,0.001,1.0,0.06);//mu2
    profile.Set(5,0.001,1.0,0.04);//mu3

    profile.Set(6,0.98,1,1.0); // mu-decay
    profile.Set(7,0.0,1.0,0.8);   // pow-decay
    profile.Set(8,0.0,1.0,0.8);   // pow-decay
    profile.Set(9,0.0,1.0,0.8);   // pow-decay
    profile.Set(10,0.0001,0.008,0.002);//mu-mix
    profile.Set(11,0.8,0.9999,0.95);//mu-mix-beta

    profile.Set(12,0.99,0.9999,0.998);
    profile.Set(13,0.001,10.0,0.001);

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

    profile.Set(24,4,mo_lpc,16);//ols-order
    profile.Set(25,4,mo_lpc,16);//ols-order
    profile.Set(26,0,mo_lpc,8);//ols-order
    profile.Set(27,0,mo_lpc,8);//ols-order

    profile.Set(28,256,4096,1280);
    profile.Set(29,32,1024,256);
    profile.Set(30,4,256,32);

    profile.Set(31,256,4096,1280);
    profile.Set(32,32,1024,256);
    profile.Set(33,4,256,32);

    profile.Set(34,0,1,0.6);
    profile.Set(35,0.1,2,0.8);
    profile.Set(36,0,10,2);

    profile.Set(37,2,32,4);
    profile.Set(38,2,32,4);

    profile.Set(39,0.98,1,1.0); // mu-decay
    profile.Set(40,0.98,1,1.0); // mu-decay

    profile.Set(41,0,1.0,1.0); // nu-mix0
    profile.Set(42,0,1.0,1.0); // nu-mix1

    profile.Set(43,0.001,0.005,0.0015);//bc-mu0
    profile.Set(44,0.001,0.005,0.0015);//bc-mu1

    profile.Set(45,4,10,5);//bias scale in bits

    profile.Set(46,0.98,1,1.0); // mu_decay
    profile.Set(47,0.98,1,1.0); // mu_decay
    profile.Set(48,0.98,1,1.0); // mu_decay
    profile.Set(49,0.98,1,1.0); // mu_decay

    profile.Set(50,0.0,1.0,0.8); //pow_decay
    profile.Set(51,0.0,1.0,0.8); //pow_decay
    profile.Set(52,0.0,1.0,0.8); //pow_decay

}
