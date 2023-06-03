#ifndef CMDLINE_H
#define CMDLINE_H

#include <string>
#include "file/wav.h"
#include "libsac/libsac.h"

const char SACHelp[] =
"usage: sac [--options] input output\n\n"
"  --help           print this message\n"
"  --encode         encode input.wav to output.sac (default)\n"
"   --normal        normal compression (default)\n"
"   --high          high compression, slow\n"
"   --optimize=#  enable frame-based optimization\n"
"     #=fast|normal|high|veryhigh|insane\n"
"   --sparse-pcm    enable pcm modelling\n"
"  --decode         decode input.sac to output.wav\n"
"  --list           list info about input.[wav|sac]\n"
"  --listfull       verbose info about input\n\n";

class CmdLine {
  enum CMODE {ENCODE,DECODE,LIST,LISTFULL};
  public:
    CmdLine();
    int Parse(int argc,char *argv[]);
    int ReadConfig(const std::string &fname);
    int Process();
  private:
    void Split(const std::string &str,std::string &key,std::string &val,const char splitval='=');
    void PrintInfo(const AudioFile &myWav);
    std::string sinputfile,soutputfile;
    CMODE mode;
    FrameCoder::coder_ctx opt;
};


#endif // CMDLINE_H
