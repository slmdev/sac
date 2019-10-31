#include "cmdline.h"

int main(int argc,char *argv[])
{
  std::cout << "Sac v0.5.1 - Lossless Audio Coder (c) Sebastian Lehmann\n";
  std::cout << "compiled on " << __DATE__ << " ";
  #ifdef __x86_64
    std::cout << "(64-bit)";
  #else
    std::cout << "(32-bit)";
  #endif
  std::cout << " with GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
  std::cout << "\n\n";
  /*int xval=32751;
  for (int bpn=0;bpn<16;bpn++) {
    std::cout << '<' << ((xval>>bpn)<<bpn) << ' ' << (xval & (~((1<<bpn)-1))) << ">\n";
  }
  return 0;*/
  CmdLine cmdline;
  int error=cmdline.Parse(argc,argv);
  if (error==0) error=cmdline.Process();
  return error;
}
