#include "cmdline.h"

int main(int argc,char *argv[])
{
  std::cout << "Sac v0.6.3 - Lossless Audio Coder (c) Sebastian Lehmann\n";
  std::cout << "compiled on " << __DATE__ << " ";
  #ifdef __x86_64
    std::cout << "(64-bit)";
  #else
    std::cout << "(32-bit)";
  #endif
  #ifdef __clang__
    std::cout << " with clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
  #elif __GNUC__ // __clang__
    std::cout << " with gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
  #endif
  std::cout << "\n\n";
  CmdLine cmdline;
  int error=cmdline.Parse(argc,argv);
  if (error==0) error=cmdline.Process();
  return error;
}
