#include "cmdline.h"
#include "info.h"

#include <iostream>

void SACInfo() {
  const int kValueWidth = 35; // 右对齐字段宽度（根据需求调整）

  std::cout << "+----------------------- Sac -----------------------+\n"
            << "| State-of-the-art lossless audio compression model |\n"
            << "+---------------------------------------------------+\n"
            << "| Copyright (c) 2024 Sebastian Lehmann  MIT License |\n"
            << "+---------------------------------------------------+\n"
            << "| Compiler      " << std::setw(kValueWidth) << std::right
            << COMPILER << " |\n"
            << "| Architecture  " << std::setw(kValueWidth) << std::right
            << ARCHITECTURE << " |\n"
            << "| AVX2 State    " << std::setw(kValueWidth) << std::right
            << AVX2_STATE << " |\n"
            << "+---------------------------------------------------+\n";
}

int main(int argc, char* argv[]) {
  SACInfo();
  CmdLine cmdline;
  int error = cmdline.Parse(argc, argv);
  if(error == 0) error = cmdline.Process();
  return error;
}