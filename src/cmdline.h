#ifndef CMDLINE_H
#define CMDLINE_H

#include <string>
#include "file/wav.h"

const char SACHelp[] =
"usage: sac [--options] input output\n\n"
"  --encode         encode input.wav to output.sac (default)\n"
"   --fast          fast, sacrifice compression\n"
"   --normal        normal compression (default)\n"
"   --high          high compression, slow\n"
"  --decode         decode input.sac to output.wav\n"
"  --list           list info about input.[wav|sac]\n"
"  --listfull       verbose info about input\n\n"
"modify encoding presets:\n"
"  --sparse-pcm     exploit sparsity of pcm-spectrum\n";

class CmdLine {
  enum CMODE {ENCODE,DECODE,LIST,LISTFULL};
  public:
    CmdLine();
    int Parse(int argc,char *argv[]);
    int Process();
  private:
    void PrintInfo(const AudioFile &myWav);
    std::string sinputfile,soutputfile;
    CMODE mode;
    int profile;
};


#endif // CMDLINE_H
