#ifndef CMDLINE_H
#define CMDLINE_H

#include <string>
#include "file/wav.h"
#include "libsac/libsac.h"

const char SACHelp[] =
"usage: sac [--options] input output\n\n"
"  --encode         encode input.wav to output.sac (default)\n"
"    --normal       normal compression (default)\n"
"    --high         high compression, slow\n"
"    --veryhigh     very high compression, really slow\n"
"    --best         you asked for it\n"
"  --decode         decode input.sac to output.wav\n"
"  --list           list info about input.sac\n"
"  --listfull       verbose info about input\n\n"
"  supported types: 1-16 bit, mono/stereo PCM\n"
"  advanced options (automatically set)\n"
"   --reset-opt     reset opt params at frame boundaries\n"
"   --optimize=#    frame-based optimization\n"
"     #=fast|normal|high|veryhigh|insane\n"
"   --sparse-pcm    enable PCM modelling\n";

class CmdLine {
  enum CMODE {ENCODE,DECODE,LIST,LISTFULL};
  public:
    CmdLine();
    int Parse(int argc,char *argv[]);
    int ReadConfig(const std::string &fname);
    int Process();
  private:
    void PrintMode();
    void PrintWav(const AudioFile &myWav);
    void Split(const std::string &str,std::string &key,std::string &val,const char splitval='=');
    std::string sinputfile,soutputfile;
    CMODE mode;
    FrameCoder::coder_ctx opt;
};


#endif // CMDLINE_H
