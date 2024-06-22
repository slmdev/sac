#include "cmdline.h"
#include "common/utils.h"
#include "common/timer.h"
#include "file/sac.h"
#include <cstring>

CmdLine::CmdLine()
:mode(ENCODE)
{
  opt.profile=0;
  opt.optimize=0;
  opt.sparse_pcm=0;
  opt.reset_profile=0;

  opt.dds_search_radius=0.15;
  opt.optimize_cost=opt.SearchCost::Entropy;
  opt.optimize_search=opt.SearchMethod::DDS;
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
  std::cout << "  Profile: ";
  if (opt.profile==0) std::cout << "normal";
  else if (opt.profile==1) std::cout << "high";
  if (opt.optimize) {
      std::cout << " (" << std::format("{:.1f}%", opt.optimize_fraction*100.0);
      std::cout << ",n=" << opt.optimize_maxnfunc<<")";
  }
  if (opt.sparse_pcm) std::cout << ", sparse-pcm";
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


int CmdLine::ReadConfig(const std::string &fname)
{
  struct profile {
    std::vector <float> ols_lambda;
    std::vector <float> ols_nu;
    std::vector <float> lms_mu0;
    std::vector <float> lms_mu1;
    std::vector <float> lms_mu2;
    std::vector <float> lms_mu_decay;
    std::vector <float> lms_pow_decay;
    std::vector <float> lms_mu_mix;
  } myprofile;
  const std::string whites=" \t\n";
  std::ifstream ifile(fname);

  if (ifile.is_open()) {
     std::string line,s;
     while (std::getline(ifile,line)) {
       s=line;
       auto firstPos = s.find_first_of("#"); // remove comment
       if (firstPos!=std::string::npos) {
         s=s.substr(0,firstPos);
       }
       StrUtils::RemoveWhite(s,whites);
       if (s.length()) {

          StrUtils::StrUpper(s);
          auto splitPos=s.find_first_of("=");
          if (splitPos!=std::string::npos) {
            //process key=val pair
            std::string key,val;
            key=s.substr(0,splitPos);
            val=s.substr(splitPos+1);
            StrUtils::RemoveWhite(val,whites);
            StrUtils::RemoveWhite(key,whites);
            if (key=="OLS_LAMBDA") StrUtils::SplitFloat(val,myprofile.ols_lambda);
            else if (key=="OLS_NU") StrUtils::SplitFloat(val,myprofile.ols_nu);
            else if (key=="LMS_MU0") StrUtils::SplitFloat(val,myprofile.lms_mu0);
            else if (key=="LMS_MU1") StrUtils::SplitFloat(val,myprofile.lms_mu1);
            else if (key=="LMS_MU2") StrUtils::SplitFloat(val,myprofile.lms_mu2);
            else if (key=="LMS_MU_DECAY") StrUtils::SplitFloat(val,myprofile.lms_mu_decay);
            else if (key=="LMS_POW_DECAY") StrUtils::SplitFloat(val,myprofile.lms_pow_decay);
            else if (key=="LMS_MU_MIX") StrUtils::SplitFloat(val,myprofile.lms_mu_mix);
            else std::cerr << "  warning: unknown option: '" << line << "'\n";
            //std::cout << "'" << s << "'\n";
          } else std::cerr << "  warning: unknown option: '" << line << "'\n";
       }
     }
    for (auto &x:myprofile.ols_lambda) std::cout << x << ' ';

    opt.profiledata.Init(12,0);
    if (myprofile.ols_lambda.size()) opt.profiledata.Set(0,myprofile.ols_lambda);
    if (myprofile.ols_nu.size()) opt.profiledata.Set(1,myprofile.ols_nu);

    std::cout << '\n';for (auto &x:myprofile.ols_nu) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_mu0) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_mu1) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_mu2) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_mu_decay) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_pow_decay) std::cout << x << ' ';
    std::cout << '\n';for (auto &x:myprofile.lms_mu_mix) std::cout << x << ' ';
              //std::cout << "'" << val << "' '" << key << "'\n";
    return 0;
  } else return 1;
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
    param=uparam=argv[k];
    StrUtils::StrUpper(uparam);
    std::string key,val;
    Split(uparam,key,val);

    if (param.length()>1 && (param[0]=='-' && param[1]=='-'))
    {
       if (key=="--ENCODE") mode=ENCODE;
       else if (key=="--DECODE") mode=DECODE;
       else if (key=="--LIST") mode=LIST;
       else if (key=="--LISTFULL") mode=LISTFULL;
       else if (key=="--VERBOSE") {
          if (val.length()) opt.verbose_level=stoi(val);
          else opt.verbose_level=1;
       }
       else if (key=="--NORMAL") opt.profile=0;
       else if (key=="--HIGH") {
         opt.profile=1;
         opt.optimize=1;
         opt.optimize_fraction=0.20;
         opt.optimize_maxnfunc=150;
         opt.sparse_pcm=0;
       } else if (key=="--VERYHIGH") {
         opt.profile=1;
         opt.optimize=1;
         opt.optimize_fraction=0.50;
         opt.optimize_maxnfunc=1000;
         opt.sparse_pcm=1;
       } else if (key=="--BEST") {
         opt.profile=1;
         opt.optimize=1;
         opt.optimize_maxnfunc=1500;
         opt.optimize_fraction=0.70;
         opt.sparse_pcm=1;
       } else if (key=="--RESET-OPT") {
         opt.reset_profile=1;
       } else if (key=="--OPTIMIZE") {
         if (val=="NONE") opt.optimize=0;
         else {
          std::vector<std::string>vs;
          StrUtils::SplitToken(val,vs,",");
          if (vs.size()>=2)  {
            opt.optimize_fraction=std::clamp(std::stod(vs[0]),0.,1.);
            opt.optimize_maxnfunc=std::clamp(std::stoi(vs[1]),0,10000);
            if (opt.optimize_fraction>0. && opt.optimize_maxnfunc>0) opt.optimize=1;
          } else std::cerr << "unknown option: " << val << '\n';
         }
       }
       else if (key=="--SPARSE-PCM") {
          if (val.length()) opt.sparse_pcm=stoi(val);
          else opt.sparse_pcm=1;
       } else if (key=="--STEREO-MS") {
        opt.stereo_ms=1;
       } else std::cout << "warning: unknown option '" << param << "'\n";
    } else {
       if (first) {sinputfile=param;first=false;}
       else soutputfile=param;
    }
    k++;
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

         bool fsupp=(myWav.getBitsPerSample()<=16) && ( (myWav.getNumChannels()==2) || (myWav.getNumChannels()==1));
         if (!fsupp)
         {
            std::cerr << "unsupported input format" << std::endl;
            std::cerr << "must be 1-16 bit, mono/stereo, PCM" << std::endl;
            myWav.Close();
            return 1;
         }
         Sac mySac(myWav);
         std::cout << "Create: '" << soutputfile << "': ";
         if (mySac.OpenWrite(soutputfile)==0) {
           std::cout << "ok\n";
           PrintMode();
           Codec myCodec;

           Timer time;
           time.start();
           myCodec.EncodeFile(myWav,mySac,opt);
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
      if (mySac.ReadHeader()==0) {
        uint8_t md5digest[16];
        mySac.ReadMD5(md5digest);

        double bps=(static_cast<double>(FileSizeSAC)*8.0)/static_cast<double>(mySac.getNumSamples()*mySac.getNumChannels());
        int kbps=round((mySac.getSampleRate()*mySac.getNumChannels()*bps)/1000);
        mySac.setKBPS(kbps);
        PrintWav(mySac);
        std::cout << "  Profile: ";
        switch (mySac.GetProfile()) {
          case 0:std::cout << "Normal";break;
          case 1:std::cout << "High";break;
        }
        std::cout << std::endl;
        std::cout << "  Ratio:   " << std::fixed << std::setprecision(3) << bps << " bits per sample\n\n";
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
            Codec myCodec;
            myCodec.DecodeFile(mySac,myWav);
            MD5::Finalize(&myWav.md5ctx);
            bool md5diff=std::memcmp(myWav.md5ctx.digest, md5digest, 16);
            std::cout << '\n';
            std::cout << "  Audio MD5: ";
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
