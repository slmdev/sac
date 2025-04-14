#include "cmdline.h"
#include "common/utils.h"
#include "common/timer.h"
#include "file/sac.h"
#include <cstring>

CmdLine::CmdLine()
:mode(ENCODE)
{
}

void CmdLine::PrintWav(const AudioFile &myWav)
{
  std::cout << "  WAVE  Codec: PCM (" << myWav.getKBPS() << " kbps)\n";
  std::cout << "  " << myWav.getSampleRate() << "Hz " << myWav.getBitsPerSample() << " Bit  ";
  if (myWav.getNumChannels()==1) std::cout << "Mono";
  else if (myWav.getNumChannels()==2) std::cout << "Stereo";
  else std::cout << myWav.getNumChannels() << " Channels";
  std::cout << "\n";
  std::cout << "  " << myWav.getNumSamples() << " Samples [" << miscUtils::getTimeStrFromSamples(myWav.getNumSamples(),myWav.getSampleRate()) << "]\n";
}

void CmdLine::PrintMode()
{
  const FrameCoder::toptim_cfg &ocfg = opt.ocfg;
  std::cout << "  Profile: ";
  std::cout << "mt" << opt.mt_mode;
  std::cout << " " << opt.max_framelen << "s";
  if (opt.adapt_block) std::cout << " ab";
  if (opt.zero_mean) std::cout << " zero-mean";
  if (opt.sparse_pcm) std::cout << " sparse-pcm";
  std::cout << '\n';
  if (opt.optimize) {
      std::cout << "  Optimize: " << std::format("{:.1f}%", ocfg.fraction*100.0);
      std::cout << ",n=" << ocfg.maxnfunc;
      std::string cost_str;
      switch (ocfg.optimize_cost) {
        case FrameCoder::SearchCost::L1:cost_str="L1";break;
        case FrameCoder::SearchCost::RMS:cost_str="rms";break;
        case FrameCoder::SearchCost::Golomb:cost_str="glb";break;
        case FrameCoder::SearchCost::Entropy:cost_str="ent";break;
        case FrameCoder::SearchCost::Bitplane:cost_str="bpn";break;
        default:break;
      }
      std::cout << "," << cost_str;
      std::cout << ",k=" << ocfg.optk;
      std::cout << '\n';
  }
  std::cout << std::endl;
}



void CmdLine::Split(const std::string &str,std::string &key,std::string &val,const char splitval)
{
  key=val="";
  std::size_t found=str.find(splitval);
  if (found!=std::string::npos) {
    key=str.substr(0,found);
    val=str.substr(found+1);
  } else {
    key=str;
  }
}

double CmdLine::stod_safe(const std::string& str)
{
    double d;
    try {
        d = std::stod(str);
    } catch (const std::invalid_argument&) {
        std::cerr << "stod: argument is invalid\n";
        throw;
    } catch (const std::out_of_range&) {
        std::cerr << "stod: argument is out of range for a double\n";
        throw;
    }
    return d;
}

int CmdLine::Parse(int argc,char *argv[])
{
  if (argc < 2) {
    std::cout << SACHelp;
    return 1;
  }

  bool first=true;
  std::string param,uparam;
  int k=1;
  while (k<argc) {
    param=argv[k];
    uparam=StrUtils::str_up(param);
    std::string key,val;
    Split(uparam,key,val);

    if (param.length()>1 && (param[0]=='-' && param[1]=='-'))
    {
       if (key=="--ENCODE") mode=ENCODE;
       else if (key=="--DECODE") mode=DECODE;
       else if (key=="--LIST") mode=LIST;
       else if (key=="--LISTFULL") mode=LISTFULL;
       else if (key=="--VERBOSE") {
          if (val.length()) opt.verbose_level=std::max(0,stoi(val));
          else opt.verbose_level=1;
       }
       else if (key=="--NORMAL") {
         opt.optimize=0;
       } else if (key=="--HIGH") {
         opt.optimize=1;
         opt.ocfg.fraction=0.075;
         opt.ocfg.maxnfunc=100;
         opt.ocfg.sigma=0.2;
         opt.ocfg.dds_cfg.c_fail_max=30;
       } else if (key=="--VERYHIGH") {
         opt.optimize=1;
         opt.ocfg.fraction=0.2;
         opt.ocfg.maxnfunc=250;
         opt.ocfg.sigma=0.2;
       } else if (key=="--EXTRAHIGH") {
         opt.optimize=1;
         opt.ocfg.fraction=0.25;
         opt.ocfg.maxnfunc=500;
         opt.ocfg.sigma=0.2;
       } else if (key=="--BEST") {
         opt.optimize=1;
         opt.ocfg.fraction=0.50;
         opt.ocfg.maxnfunc=1000;
         opt.ocfg.sigma=0.20;
         opt.ocfg.optimize_cost=FrameCoder::SearchCost::Bitplane;
       } else if (key=="--INSANE") {
         opt.optimize=1;
         opt.ocfg.fraction=0.75;
         opt.ocfg.maxnfunc=1500;
         opt.ocfg.sigma=0.25;
         opt.ocfg.optimize_cost=FrameCoder::SearchCost::Bitplane;
       } else if (key=="--OPTIMIZE") {
         if (val=="NO" || val=="0") opt.optimize=0;
         else {
          std::vector<std::string>vs;
          StrUtils::SplitToken(val,vs,",");
          if (vs.size()>=2)  {
            opt.ocfg.fraction=std::clamp(stod_safe(vs[0]),0.,1.);
            opt.ocfg.maxnfunc=std::clamp(std::stoi(vs[1]),0,50000);
            if (vs.size()>=3) {
              std::string cf=StrUtils::str_up(vs[2]);
              if (cf=="L1") opt.ocfg.optimize_cost = FrameCoder::SearchCost::L1;
              else if (cf=="RMS") opt.ocfg.optimize_cost = FrameCoder::SearchCost::RMS;
              else if (cf=="GLB") opt.ocfg.optimize_cost = FrameCoder::SearchCost::Golomb;
              else if (cf=="ENT") opt.ocfg.optimize_cost = FrameCoder::SearchCost::Entropy; //default
              else if (cf=="BPN") opt.ocfg.optimize_cost = FrameCoder::SearchCost::Bitplane;
              else std::cerr << "warning: unknown cost function '" << vs[2] << "'\n";
            }
            if (vs.size()>=4) {
              opt.ocfg.optk=std::clamp(stoi(vs[3]),1,32);
            }
            if (opt.ocfg.fraction>0. && opt.ocfg.maxnfunc>0) opt.optimize=1;
            else opt.optimize=0;
          } else std::cerr << "unknown option: " << val << '\n';
         }
       }
       else if (key=="--FRAMELEN") {
         if (val.length()) opt.max_framelen=std::max(0,stoi(val));
       }
       else if (key=="--MT-MODE") {
         if (val.length()) opt.mt_mode=std::max(0,stoi(val));
       }
       else if (key=="--SPARSE-PCM") {
          if (val=="NO" || val=="0") opt.sparse_pcm=0;
          else opt.sparse_pcm=1;
       } else if (key=="--STEREO-MS") {
         opt.stereo_ms=1;
       } else if (key=="--OPT-RESET") {
         opt.ocfg.reset=1;
       } else if (key=="--OPT-CFG") {
         std::vector<std::string>vs;
         StrUtils::SplitToken(val,vs,",");
         if (vs.size()>=1) {
            std::string val=StrUtils::str_up(vs[0]);
            if (val=="DDS") opt.ocfg.optimize_search=FrameCoder::SearchMethod::DDS;
            else if (val=="DE") opt.ocfg.optimize_search=FrameCoder::SearchMethod::DE;
            else std::cerr << "  warning: invalid val='"<<val<<"'\n";
         }
         if (vs.size()>=2) opt.ocfg.num_threads = std::clamp(std::stoi(vs[1]),0,256);
         if (vs.size()>=3) opt.ocfg.sigma=std::clamp(stod_safe(vs[2]),0.,1.);
       } else if (key=="--ADAPT-BLOCK") {
         if (val=="NO" || val=="0") opt.adapt_block=0;
         else opt.adapt_block=1;
       } else if (key=="--ZERO-MEAN") {
         if (val=="NO" || val=="0") opt.zero_mean=0;
         else opt.zero_mean=1;
       }
       else std::cerr << "warning: unknown option '" << param << "'\n";
    } else {
       if (first) {sinputfile=param;first=false;}
       else soutputfile=param;
    }
    k++;
  }
  // configure opt method
  if (opt.ocfg.optimize_search==FrameCoder::SearchMethod::DDS)
  {
    opt.ocfg.dds_cfg.nfunc_max=opt.ocfg.maxnfunc;
    opt.ocfg.dds_cfg.num_threads=opt.ocfg.num_threads; // also accepts zero
    opt.ocfg.dds_cfg.sigma_init=opt.ocfg.sigma;
  } else if (opt.ocfg.optimize_search==FrameCoder::SearchMethod::DE)
  {
    opt.ocfg.de_cfg.nfunc_max=opt.ocfg.maxnfunc;
    opt.ocfg.de_cfg.num_threads=std::max(opt.ocfg.num_threads,1);
    opt.ocfg.de_cfg.sigma_init=opt.ocfg.sigma;
  }

  return 0;
}

int CmdLine::Process()
{
  Timer myTimer;
  myTimer.start();

  if (mode==ENCODE) {
    Wav myWav(opt.verbose_level>0);
    std::cout << "Open: '" << sinputfile << "': ";
    if (myWav.OpenRead(sinputfile)==0) {
      std::cout << "ok (" << myWav.getFileSize() << " Bytes)\n";
      if (myWav.ReadHeader()==0) {
         PrintWav(myWav);

         bool fsupp=(myWav.getBitsPerSample()<=24) && ( (myWav.getNumChannels()==2) || (myWav.getNumChannels()==1));
         if (!fsupp)
         {
            std::cerr << "unsupported input format" << std::endl;
            std::cerr << "must be 1-16 bit, mono/stereo, pcm" << std::endl;
            myWav.Close();
            return 1;
         }
         Sac mySac(myWav);
         std::cout << "Create: '" << soutputfile << "': ";
         if (mySac.OpenWrite(soutputfile)==0) {
           std::cout << "ok\n";
           PrintMode();
           Codec myCodec(opt);

           Timer time;
           time.start();
           myCodec.EncodeFile(myWav,mySac);
           time.stop();

           uint64_t infilesize=myWav.getFileSize();
           uint64_t outfilesize=mySac.readFileSize();
           double r=0.,bps=0.;
           if (outfilesize) {
             r=outfilesize*100.0/infilesize;
             bps=(outfilesize*8.)/static_cast<double>(myWav.getNumSamples()*myWav.getNumChannels());
           }
           double xrate=0.0;
           if (time.elapsedS() > 0.0)
            xrate=(myWav.getNumSamples()/double(myWav.getSampleRate()))/time.elapsedS();

           std::cout << "\n  " << infilesize << "->" << outfilesize<< "=";
           std::cout << std::format("{:.1f}",r) << "% (" << std::format("{:.3f}",bps) <<" bps)";
           std::cout << "  " << std::format("{:.3f}x",xrate) << '\n';
           mySac.Close();
         } else std::cout << "could not create\n";
      } else std::cout << "warning: input is not a valid .wav file\n";
      myWav.Close();
    } else std::cout << "could not open\n";
  } else if (mode==LIST || mode==LISTFULL || mode==DECODE) {
    Sac mySac;
    std::cout << "Open: '" << sinputfile << "': ";
    if (mySac.OpenRead(sinputfile)==0) {
      std::streampos FileSizeSAC = mySac.getFileSize();
      std::cout << "ok (" << FileSizeSAC << " Bytes)\n";
      if (mySac.ReadSACHeader()==0) {
        uint8_t md5digest[16];
        mySac.ReadMD5(md5digest);
        double bps=(static_cast<double>(FileSizeSAC)*8.0)/static_cast<double>(mySac.getNumSamples()*mySac.getNumChannels());
        int kbps=round((mySac.getSampleRate()*mySac.getNumChannels()*bps)/1000);
        mySac.setKBPS(kbps);
        PrintWav(mySac);
        std::cout << "  Profile: ";
        std::cout << "mt" << opt.mt_mode;
        std::cout << " " << static_cast<int>(mySac.mcfg.max_framelen) << "s";
        std::cout << std::endl;
        std::cout << "  Ratio:   " << std::fixed << std::setprecision(3) << bps << " bps\n\n";
        std::cout << "  Audio MD5: ";
        for (auto x : md5digest) std::cout << std::hex << (int)x;
        std::cout << std::dec << '\n';


        if (mode==LISTFULL) {
          Codec myCodec;
          myCodec.ScanFrames(mySac);
        } else if (mode==DECODE) {

          Wav myWav(mySac);
          std::cout << "Create: '" << soutputfile << "': ";
          if (myWav.OpenWrite(soutputfile)==0) {
            std::cout << "ok\n";

            Timer time;
            time.start();

            Codec myCodec(opt);
            myCodec.DecodeFile(mySac,myWav);
            MD5::Finalize(&myWav.md5ctx);
            time.stop();

            double xrate=0.0;
            if (time.elapsedS() > 0.0)
            xrate=(myWav.getNumSamples()/double(myWav.getSampleRate()))/time.elapsedS();
            std::cout << "\n  Speed " << std::format("{:.3f}x",xrate) << '\n';

            std::cout << "  Audio MD5: ";
            bool md5diff=std::memcmp(myWav.md5ctx.digest, md5digest, 16);
            if (!md5diff) std::cout << "ok\n";
            else {
              std::cout << "Error (";
              for (auto x : myWav.md5ctx.digest) std::cout << std::hex << (int)x;
              std::cout << std::dec << ")\n";
            }
            myWav.Close();
          } else std::cout << "could not create\n";
        }
      } else std::cout << "warning: input is not a valid .sac file\n";
    }
  } else std::cout << "could not open\n";

  myTimer.stop();
  std::cout << "\n  Time:    [" << miscUtils::getTimeStrFromSeconds(round(myTimer.elapsedS())) << "]" << std::endl;
  return 0;
}
