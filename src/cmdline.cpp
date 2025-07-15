#include <format>
#include <cstring>
#include "cmdline.h"
#include "common/utils.h"
#include "common/timer.h"
#include "file/sac.h"


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

std::string CmdLine::CostStr(const FrameCoder::SearchCost cost_func)
{
  using enum FrameCoder::SearchCost;
  std::string rstr;
  switch (cost_func) {
    case L1:rstr="L1";break;
    case RMS:rstr="rms";break;
    case Golomb:rstr="glb";break;
    case Entropy:rstr="ent";break;
    case Bitplane:rstr="bpn";break;
    default:break;
  }
  return rstr;
}

std::string CmdLine::SearchStr(const FrameCoder::SearchMethod search_func)
{
  using enum FrameCoder::SearchMethod;
  std::string rstr;
  switch (search_func) {
    case DDS:rstr="DDS";break;
    case DE:rstr="DE";break;
    default: break;
  }
  return rstr;
}

void CmdLine::PrintMode()
{
  const FrameCoder::toptim_cfg &ocfg = cfg.ocfg;
  std::cout << "  Profile: ";
  std::cout << "mt" << cfg.mt_mode;
  std::cout << " " << cfg.max_framelen << "s";
  if (cfg.adapt_block) std::cout << " ab";
  if (cfg.zero_mean) std::cout << " zero-mean";
  if (cfg.sparse_pcm) std::cout << " sparse-pcm";
  std::cout << '\n';
  if (cfg.optimize) {
      std::cout << "  Optimize: ";
      std::cout << SearchStr(ocfg.optimize_search);
      std::cout << " " << std::format("{:.1f}%", ocfg.fraction*100.0);
      std::cout << ",n=" << ocfg.maxnfunc;
      std::cout << "," << CostStr(ocfg.optimize_cost);
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

int CmdLine::Parse(int argc,const char *argv[])
{
  if (argc < 2) {
    std::cout << SACHelp;
    return 1;
  }

  bool first=true;
  std::string param;
  int k=1;
  while (k<argc) {
    param=argv[k];
    std::string uparam=StrUtils::str_up(param);
    std::string key,val;
    Split(uparam,key,val);

    if (param.length()>1 && (param[0]=='-' && param[1]=='-'))
    {
       if (key=="--ENCODE") mode=ENCODE;
       else if (key=="--DECODE") mode=DECODE;
       else if (key=="--LIST") mode=LIST;
       else if (key=="--LISTFULL") mode=LISTFULL;
       else if (key=="--VERBOSE") {
          if (val.length()) cfg.verbose_level=std::max(0,stoi(val));
          else cfg.verbose_level=1;
       }
       else if (key=="--NORMAL") {
         cfg.optimize=0;
       } else if (key=="--HIGH") {
         cfg.optimize=1;
         cfg.ocfg.fraction=0.075;
         cfg.ocfg.maxnfunc=100;
         cfg.ocfg.sigma=0.2;
         cfg.ocfg.dds_cfg.c_fail_max=30;
       } else if (key=="--VERYHIGH") {
         cfg.optimize=1;
         cfg.ocfg.fraction=0.2;
         cfg.ocfg.maxnfunc=250;
         cfg.ocfg.sigma=0.25;
       } else if (key=="--EXTRAHIGH") {
         cfg.optimize=1;
         cfg.ocfg.fraction=0.25;
         cfg.ocfg.maxnfunc=500;
         cfg.ocfg.sigma=0.25;
       } else if (key=="--BEST") {
         cfg.optimize=1;
         cfg.ocfg.fraction=0.50;
         cfg.ocfg.maxnfunc=1000;
         cfg.ocfg.sigma=0.25;
         cfg.ocfg.optimize_cost=FrameCoder::SearchCost::Bitplane;
       } else if (key=="--INSANE") {
         cfg.optimize=1;
         cfg.ocfg.fraction=0.5;
         cfg.ocfg.maxnfunc=1500;
         cfg.ocfg.sigma=0.25;
         cfg.ocfg.optimize_cost=FrameCoder::SearchCost::Bitplane;
       } else if (key=="--OPTIMIZE") {
         if (val=="NO" || val=="0") cfg.optimize=0;
         else {
          std::vector<std::string>vs;
          StrUtils::SplitToken(val,vs,",");
          if (vs.size()>=2)  {
            cfg.ocfg.fraction=std::clamp(stod_safe(vs[0]),0.,1.);
            cfg.ocfg.maxnfunc=std::clamp(std::stoi(vs[1]),0,50000);
            if (vs.size()>=3) {
              std::string cf=StrUtils::str_up(vs[2]);
              if (cf=="L1") cfg.ocfg.optimize_cost = FrameCoder::SearchCost::L1;
              else if (cf=="RMS") cfg.ocfg.optimize_cost = FrameCoder::SearchCost::RMS;
              else if (cf=="GLB") cfg.ocfg.optimize_cost = FrameCoder::SearchCost::Golomb;
              else if (cf=="ENT") cfg.ocfg.optimize_cost = FrameCoder::SearchCost::Entropy; //default
              else if (cf=="BPN") cfg.ocfg.optimize_cost = FrameCoder::SearchCost::Bitplane;
              else std::cerr << "warning: unknown cost function '" << vs[2] << "'\n";
            }
            if (vs.size()>=4) {
              cfg.ocfg.optk=std::clamp(stoi(vs[3]),1,32);
            }
            if (cfg.ocfg.fraction>0. && cfg.ocfg.maxnfunc>0) cfg.optimize=1;
            else cfg.optimize=0;
          } else std::cerr << "unknown option: " << val << '\n';
         }
       }
       else if (key=="--FRAMELEN") {
         if (val.length()) cfg.max_framelen=std::max(0,stoi(val));
       }
       else if (key=="--MT-MODE") {
         if (val.length()) cfg.mt_mode=std::max(0,stoi(val));
       }
       else if (key=="--SPARSE-PCM") {
          if (val=="NO" || val=="0") cfg.sparse_pcm=0;
          else cfg.sparse_pcm=1;
       } else if (key=="--STEREO-MS") {
         cfg.stereo_ms=1;
       } else if (key=="--OPT-RESET") {
         cfg.ocfg.reset=1;
       } else if (key=="--OPT-CFG") {
         std::vector<std::string>vs;
         StrUtils::SplitToken(val,vs,",");
         if (vs.size()>=1) {
            std::string cval=StrUtils::str_up(vs[0]);
            if (cval=="DDS") cfg.ocfg.optimize_search=FrameCoder::SearchMethod::DDS;
            else if (cval=="DE") cfg.ocfg.optimize_search=FrameCoder::SearchMethod::DE;
            else if (cval=="CMA") cfg.ocfg.optimize_search=FrameCoder::SearchMethod::CMA;
            else std::cerr << "  warning: invalid opt='"<<cval<<"'\n";
         }
         if (vs.size()>=2) cfg.ocfg.num_threads = std::clamp(std::stoi(vs[1]),0,256);
         if (vs.size()>=3) cfg.ocfg.sigma=std::clamp(stod_safe(vs[2]),0.,1.);
       } else if (key=="--ADAPT-BLOCK") {
         if (val=="NO" || val=="0") cfg.adapt_block=0;
         else cfg.adapt_block=1;
       } else if (key=="--ZERO-MEAN") {
         if (val=="NO" || val=="0") cfg.zero_mean=0;
         else cfg.zero_mean=1;
       }
       else std::cerr << "warning: unknown option '" << param << "'\n";
    } else {
       if (first) {sinputfile=param;first=false;}
       else soutputfile=param;
    }
    k++;
  }
  // configure opt method
  if (cfg.ocfg.optimize_search==FrameCoder::SearchMethod::DDS)
  {
    cfg.ocfg.dds_cfg.nfunc_max=cfg.ocfg.maxnfunc;
    cfg.ocfg.dds_cfg.num_threads=cfg.ocfg.num_threads; // also accepts zero
    cfg.ocfg.dds_cfg.sigma_init=cfg.ocfg.sigma;
  } else if (cfg.ocfg.optimize_search==FrameCoder::SearchMethod::DE)
  {
    cfg.ocfg.de_cfg.nfunc_max=cfg.ocfg.maxnfunc;
    cfg.ocfg.de_cfg.num_threads=std::max(cfg.ocfg.num_threads,1);
    cfg.ocfg.de_cfg.sigma_init=cfg.ocfg.sigma;
  } else if (cfg.ocfg.optimize_search==FrameCoder::SearchMethod::CMA)
  {
    cfg.ocfg.cma_cfg.nfunc_max=cfg.ocfg.maxnfunc;
    cfg.ocfg.cma_cfg.num_threads=std::max(cfg.ocfg.num_threads,1);
    cfg.ocfg.cma_cfg.sigma_init=cfg.ocfg.sigma;
  }

  return 0;
}

int CmdLine::Process()
{
  Timer myTimer;
  myTimer.start();

  if (mode==ENCODE) {
    Wav myWav(cfg.verbose_level>0);
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
           Codec myCodec(cfg);

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
        std::cout << "mt" << cfg.mt_mode;
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

            Codec myCodec(cfg);
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
