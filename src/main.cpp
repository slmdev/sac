#include "cmdline.h"

/*class DDS {
  public:
    DDS(int maxk)
    :maxk(maxk)
    {

    };
  private:
    int maxk;
};




class Func2D {
  public:
    Func2D(){hist[0]=hist[1]=1;};
    double Next() {
      double t=exp(-hist[0]*hist[0]);
      double r=(0.8-0.5*t)*hist[0]-
      (0.3+0.9*t)*hist[1]+0.1*sin(3.1415*hist[0]);

      hist[1]=hist[0];hist[0]=r;
      return r;
    }
  double hist[2];
};*/


int main(int argc,char *argv[])
{
  std::cout << "Sac v0.0.7a - Lossless Audio Coder (c) Sebastian Lehmann\n";
  std::cout << "compiled on " << __DATE__ << " ";
  #ifdef __x86_64
    std::cout << "(64-bit)";
  #else
    std::cout << "(32-bit)";
  #endif
  std::cout << "\n\n";

  int error=0;
  CmdLine cmdline;
  error=cmdline.Parse(argc,argv);
  if (error==0) error=cmdline.Process();
  return error;
}
