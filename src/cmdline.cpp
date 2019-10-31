#include "cmdline.h"
#include "common/utils.h"
#include "common/timer.h"
#include "file/sac.h"

CmdLine::CmdLine()
:mode(ENCODE)
{
  opt.profile=0;
  opt.optimize=0;
  opt.sparse_pcm=0;
  opt.optimize_mode=0;
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
  //std::cout << "'" << key << "':'" << val << "'\n";
  //system("Pause");
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
    if (param.length()>1 && (param[0]=='-' && param[1]=='-')) {
       if (key=="--ENCODE") mode=ENCODE;
       else if (key=="--DECODE") mode=DECODE;
       else if (key=="--LIST") mode=LIST;
       else if (key=="--LISTFULL") mode=LISTFULL;
       else if (key=="--NORMAL") opt.profile=0;
       else if (key=="--HIGH") opt.profile=1;
       else if (key=="--OPTIMIZE") {
         opt.optimize=1;
         if (val=="FAST") {
           opt.optimize_mode=0;
         } else if (val=="NORMAL") {
           opt.optimize_mode=1;
         } else if (val=="HIGH") {
           opt.optimize_mode=2;
         } else if (val=="INSANE") {
           opt.optimize_mode=3;
         }
       }
       else if (key=="--SPARSE-PCM") opt.sparse_pcm=1;
       else std::cout << "warning: unknown option '" << param << "'\n";
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

  /*std::string fname_cfg="sac.cfg";
  if (ReadConfig(fname_cfg)>0) {
    std::cerr << "could not open: " << fname_cfg << '\n';
    return 1;
  } else {
    return 0;
  }*/

  if (mode==ENCODE) {
    Wav myWav(false);
    std::cout << "Open: '" << sinputfile << "': ";
    if (myWav.OpenRead(sinputfile)==0) {
      std::cout << "ok (" << myWav.getFileSize() << " Bytes)\n";
      if (myWav.ReadHeader()==0) {
         if (myWav.getBitsPerSample()!=16 || myWav.getNumChannels()!=2) {
            std::cout << "Unsupported input format." << std::endl;
            myWav.Close();
            return 1;
         }
         PrintInfo(myWav);
         Sac mySac(myWav);
         std::cout << "Create: '" << soutputfile << "': ";
         if (mySac.OpenWrite(soutputfile)==0) {
           std::cout << "ok\n";
           std::cout << "  Profile: ";
           if (opt.profile==0) std::cout << "Normal";
           else if (opt.profile==1) std::cout << "High";
           if (opt.optimize) {
            switch (opt.optimize_mode) {
              case 0: std::cout << " (optimize fast)";break;
              case 1: std::cout << " (optimize normal)";break;
              case 2: std::cout << " (optimize high)";break;
              case 3: std::cout << " (optimize insane)";break;
              default: std::cout << " (unknown profile)";break;
            }
           }
           std::cout << std::endl;
           Codec myCodec;

           myCodec.EncodeFile(myWav,mySac,opt);

           uint64_t infilesize=myWav.getFileSize();
           uint64_t outfilesize=mySac.readFileSize();
           double r=0.,bps=0.;
           if (outfilesize) {
             r=outfilesize*100.0/infilesize;
             bps=(outfilesize*8.)/static_cast<double>(myWav.getNumSamples()*myWav.getNumChannels());
           }
           std::cout << std::endl << "  " << infilesize << "->" << outfilesize<< "=" <<miscUtils::ConvertFixed(r,1) << "% (" << miscUtils::ConvertFixed(bps,3)<<" bps)"<<std::endl;
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
        double bps=(static_cast<double>(FileSizeSAC)*8.0)/static_cast<double>(mySac.getNumSamples()*mySac.getNumChannels());
        int kbps=round((mySac.getSampleRate()*mySac.getNumChannels()*bps)/1000);
        mySac.setKBPS(kbps);
        PrintInfo(mySac);
        std::cout << "  Profile: ";
        switch (mySac.GetProfile()) {
          case 0:std::cout << "Normal";break;
          case 1:std::cout << "High";break;
        }
        std::cout << std::endl;
        std::cout << "  Ratio:   " << std::fixed << std::setprecision(3) << bps << " bits per sample\n\n";


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
            myWav.Close();
          } else std::cout << "could not create\n";
        }
      } else std::cout << "warning: input is not a valid .sac file\n";
    }
  } else std::cout << "could not open\n";

  myTimer.stop();
  std::cout << "\nTotal time: [" << miscUtils::getTimeStrFromSeconds(round(myTimer.elapsedS())) << "]" << std::endl;
  return 0;
}
