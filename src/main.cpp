#include "cmdline.h"
#include <variant>

#define SAC_VERSION "0.7.11"

#include <thread>
#include <future>
int test(int ch)
{
  return ch*2;
}

int main(int argc,char *argv[])
{
  std::cout << "Sac v" << SAC_VERSION << " - Lossless Audio Coder (c) Sebastian Lehmann\n";
  std::cout << "compiled on " << __DATE__ << " ";
  #ifdef __x86_64
    std::cout << "(64-bit)";
  #else
    std::cout << "(32-bit)";
  #endif
  #ifdef __clang__
    std::cout << " clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
  #elif __GNUC__ // __clang__
    std::cout << " gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
  #endif
  std::cout << "\n";

  CmdLine cmdline;
  int error=cmdline.Parse(argc,argv);
  if (error==0) error=cmdline.Process();
}
