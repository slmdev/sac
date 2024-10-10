#ifndef CMDLINE_H
#define CMDLINE_H

#include <string>
#include "file/wav.h"
#include "libsac/libsac.h"

const char SACHelp[] =
"usage: sac [--options] input output\n\n"
"  --encode            encode input.wav to output.sac (def)\n"
"    --normal          normal compression (def)\n"
"    --high            high compression, slow\n"
"    --veryhigh        very high compression, really slow\n"
"    --best            you asked for it\n"
"    --insane          :>\n"
"  --decode            decode input.sac to output.wav\n"
"  --list              list info about input.sac\n"
"  --listfull          verbose info about input\n"
"  --verbose=n         verbosity level (def=0)\n\n"
"  supported types: 1-16 bit, mono/stereo pcm\n"
"  advanced options    (automatically set)\n"
"   --reset-opt        reset opt params at frame boundaries\n"
"   --optimize=#       frame-based optimization\n"
"     no|s,n(,c)       s=[0,1.0],n=[0,10000]\n"
"                      c=[l1,rms,glb,ent,bpn]\n"
"   --mt-mode=n        multi-threading level n=[0-2]\n"
"   --zero-mean        zero-mean input\n"
"   --framelen=n       def=8 (seconds)\n"
"   --sparse-pcm       enable pcm modelling\n";

class CmdLine {
  enum CMODE {ENCODE,DECODE,LIST,LISTFULL};
  public:
    CmdLine();
    int Parse(int argc,char *argv[]);
    int Process();
  private:
    double stod_safe(const std::string& str);
    void PrintMode();
    void PrintWav(const AudioFile &myWav);
    void Split(const std::string &str,std::string &key,std::string &val,const char splitval='=');
    std::string sinputfile,soutputfile;
    CMODE mode;
    FrameCoder::coder_ctx opt;
};


#endif // CMDLINE_H
