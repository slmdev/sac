#include "cmdline.h"
#include "common/utils.h"
#include "common/timer.h"
#include "file/sac.h"
#include "libsac/libsac.h"

CmdLine::CmdLine()
:mode(ENCODE),profile(1)
{

}

void CmdLine::PrintInfo(const AudioFile &myWav)
{
  std::cout << "  WAVE  Codec: PCM (" << myWav.getKBPS() << " kbps)\n";
  std::cout << "  " << myWav.getSampleRate() << "Hz " << myWav.getBitsPerSample() << " Bit  ";
  if (myWav.getNumChannels()==1) std::cout << "Mono";
  else if (myWav.getNumChannels()==2) std::cout << "Stereo";
  else std::cout << myWav.getNumChannels() << " Channels";
  std::cout << "\n";
  std::cout << "  " << myWav.getNumSamples() << " Samples [" << miscUtils::getTimeStrFromSamples(myWav.getNumSamples(),myWav.getSampleRate()) << "]\n";
}

int CmdLine::Parse(int argc,char *argv[])
{
  if (argc < 2) {
    std::cout << SACHelp;
    return 1;
  }

  bool first=true;
  std::string param,uparam;
  for (int i=1;i<argc;i++) {
    param=uparam=argv[i];
    strUtils::strUpper(uparam);
    if (param.length()>1 && (param[0]=='-' && param[1]=='-')) {
       if (uparam.compare("--ENCODE")==0) mode=ENCODE;
       else if (uparam.compare("--DECODE")==0) mode=DECODE;
       else if (uparam.compare("--LIST")==0) mode=LIST;
       else if (uparam.compare("--LISTFULL")==0) mode=LISTFULL;
       else if (uparam.compare("--FAST")==0) profile=0;
       else if (uparam.compare("--NORMAL")==0) profile=1;
       else if (uparam.compare("--HIGH")==0) profile=2;
       else std::cout << "warning: unknown option '" << param << "'\n";
    } else {
       if (first) {sinputfile=param;first=false;}
       else soutputfile=param;
    }
  }
  return 0;
}

int CmdLine::Process()
{
  Timer myTimer;
  myTimer.start();

  if (mode==ENCODE) {
    Wav myWav(true);
    std::cout << "Open: '" << sinputfile << "': ";
    if (myWav.OpenRead(sinputfile)==0) {
      std::cout << "ok (" << myWav.getFileSize() << " Bytes)\n";
      if (myWav.ReadHeader()==0) {
         if (myWav.getBitsPerSample()!=16 || myWav.getNumChannels()!=2) {
            std::cout << "Unsupported input format." << std::endl;
            myWav.Close();
            return 1;
         }
         Sac mySac(myWav);
         std::cout << "Create: '" << soutputfile << "': ";
         if (mySac.OpenWrite(soutputfile)==0) {
           std::cout << "ok\n";
           PrintInfo(myWav);
           std::cout << "  Profile: ";
           if (profile==0) std::cout << "Fast" << std::endl;
           else if (profile==1) std::cout << "Normal" << std::endl;
           else if (profile==2) std::cout << "High" << std::endl;
           Codec myCodec;
           myCodec.EncodeFile(myWav,mySac,profile);
           uint64_t infilesize=myWav.getFileSize();
           uint64_t outfilesize=mySac.readFileSize();
           double r=0.,bps=0.;
           if (outfilesize) {
             r=outfilesize*100./infilesize;
             bps=(outfilesize*8.)/static_cast<double>(myWav.getNumSamples()*myWav.getNumChannels());
           }
           std::cout << std::endl << "  " << infilesize << "->" << outfilesize<< "=" <<miscUtils::ConvertFixed(r,2) << "% (" << miscUtils::ConvertFixed(bps,2)<<" bps)"<<std::endl;

           mySac.Close();
         } else std::cout << "could not create\n";
      } else std::cout << "warning: input is not a valid .wav file\n";
      myWav.Close();
    } else std::cout << "could not open\n";
  } else if (mode==DECODE) {
    std::cout << "Open: '" << sinputfile << "': ";
    Sac mySac;
    if (mySac.OpenRead(sinputfile)==0) {
      std::cout << "ok (" << mySac.getFileSize() << " Bytes)\n";
      if (mySac.ReadHeader()==0) {
           PrintInfo(mySac);
           Wav myWav(mySac,true);
           std::cout << "Create: '" << soutputfile << "': ";
           if (myWav.OpenWrite(soutputfile)==0) {
             std::cout << "ok\n";
             std::cout << "  Profile: ";
             switch (mySac.GetProfile()) {
               case 0:std::cout << "Fast";break;
               case 1:std::cout << "Normal";break;
               case 2:std::cout << "High";break;
             }
             std::cout << std::endl;
             Codec myCodec;
             myCodec.DecodeFile(mySac,myWav);
             myWav.Close();
           } else std::cout << "could not create\n";
      } else std::cout << "warning: input is not a valid .sac file\n";
      mySac.Close();
    } else std::cout << "could not open\n";
  } else if (mode==LIST || mode==LISTFULL) {
    Sac mySac;
    std::cout << "Open: '" << sinputfile << "': ";
    if (mySac.OpenRead(sinputfile)==0) {
      std::cout << "ok (" << mySac.getFileSize() << " Bytes)\n";
      if (mySac.ReadHeader()==0) {
        PrintInfo(mySac);
        std::cout << "  Profile: ";
        switch (mySac.GetProfile()) {
          case 0:std::cout << "Fast";break;
          case 1:std::cout << "Normal";break;
          case 2:std::cout << "High";break;
        }
        std::cout << std::endl;
        if (mode==LISTFULL) {
          Codec myCodec;
          myCodec.ScanFrames(mySac);
        }
      } else std::cout << "warning: input is not a valid .sac file\n";
    }
  } else std::cout << "could not open\n";

  myTimer.stop();
  std::cout << "\nTotal time: [" << miscUtils::getTimeStrFromSeconds(round(myTimer.elapsedS())) << "]" << std::endl;
  return 0;
}
