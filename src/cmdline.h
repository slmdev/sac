#ifndef CMDLINE_H
#define CMDLINE_H

#include <string>
#include "file/wav.h"
#include "libsac/libsac.h"

const char SACHelp[] =
"usage: sac [--options] input output\n\n"
"  --encode            encode input.wav to output.sac (def)\n"
"    --normal|high|veryhigh|extrahigh compression (def=normal)\n"
"    --best            you asked for it\n\n"
"  --decode            decode input.sac to output.wav\n"
"  --list              list info about input.sac\n"
"  --listfull          verbose info about input\n"
"  --verbose           verbose output\n\n"
"  supported types: 1-16 bit, mono/stereo pcm\n"
"  advanced options    (automatically set)\n"
"   --optimize=#       frame-based optimization\n"
"     no|s,n,c,k       s=[0,1.0],n=[0,10000]\n"
"                      c=[l1,rms,glb,ent,bpn] k=[1,32]\n"
"   --opt-cfg=#        configure optimization method\n"
"     de|dds,nt,s      nt=num threads,s=search radius (def=0.2)\n"
"   --opt-reset        reset opt params at frame boundaries\n"
"   --mt-mode=n        multi-threading level n=[0-2]\n"
"   --zero-mean        zero-mean input\n"
"   --adapt-block      adaptive frame splitting\n"
"   --framelen=n       def=20 seconds\n"
"   --sparse-pcm       enable pcm modelling\n";

class CmdLine {
  enum CMODE {ENCODE,DECODE,LIST,LISTFULL};
  public:
    CmdLine();
    int Parse(int argc,const char *argv[]);
    int Process();
  private:
    double stod_safe(const std::string& str);
    std::string CostStr(const FrameCoder::SearchCost cost_func);
    std::string SearchStr(const FrameCoder::SearchMethod search_func);
    void PrintMode();
    void PrintWav(const AudioFile &myWav);
    void Split(const std::string &str,std::string &key,std::string &val,const char splitval='=');
    std::string sinputfile,soutputfile;
    CMODE mode;
    FrameCoder::coder_ctx opt;
};


#endif // CMDLINE_H
