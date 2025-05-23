#ifndef CODEC_H
#define CODEC_H

#include "../file/wav.h"
#include "../file/sac.h"
#include "cost.h"
#include "profile.h"
#include "../opt/dds.h"
#include "../opt/de.h"

class FrameCoder {
  public:
    enum SearchCost {L1,RMS,Entropy,Golomb,Bitplane};
    enum SearchMethod {DDS,DE};

    typedef std::vector <std::vector<int32_t>> tch_samples;

    struct toptim_cfg {
      OptDDS::DDSCfg dds_cfg;
      OptDE::DECfg de_cfg;
      int reset=0;
      double fraction=0;
      int maxnfunc=0;
      int num_threads=0;
      double sigma=0.2;
      int optk=4;
      SearchMethod optimize_search=SearchMethod::DDS;
      SearchCost optimize_cost=SearchCost::Entropy;
    };
    struct coder_ctx {
      int optimize=0;
      int sparse_pcm=1;
      int zero_mean=1;
      int max_framelen=20;
      int verbose_level=0;
      int stereo_ms=0;
      int mt_mode=2;
      int adapt_block=1;

      toptim_cfg ocfg;
      SacProfile profiledata;
    };
    FrameCoder(int numchannels,int framesize,const coder_ctx &opt);
    void SetNumSamples(int nsamples){numsamples_=nsamples;};
    int GetNumSamples(){return numsamples_;};
    void Predict();
    void Unpredict();
    void Encode();
    void Decode();
    void WriteEncoded(AudioFile &fout);
    void ReadEncoded(AudioFile &fin);
    std::vector <std::vector<int32_t>>samples,error,s2u_error,s2u_error_map,pred;
    std::vector <BufIO> encoded,enc_temp1,enc_temp2;
    std::vector <SacProfile::FrameStats> framestats;

    static int WriteBlockHeader(std::fstream &file, const std::vector<SacProfile::FrameStats> &framestats, int ch);
    static int ReadBlockHeader(std::fstream &file, std::vector<SacProfile::FrameStats> &framestats, int ch);
  private:
    void CnvError_S2U(const tch_samples &error,int numsamples);
    void SetParam(Predictor::tparam &param,const SacProfile &profile,bool optimize=false);
    void PrintProfile(SacProfile &profile);
    void EncodeProfile(const SacProfile &profile,std::vector <uint8_t>&buf);
    void DecodeProfile(SacProfile &profile,const std::vector <uint8_t>&buf);
    void AnalyseMonoChannel(int ch, int numsamples);
    double AnalyseStereoChannel(int ch0, int ch1, int numsamples);
    void ApplyMs(int ch0, int ch1, int numsamples);
    //void InterChannel(int ch0,int ch1,int numsamples);
    int EncodeMonoFrame_Normal(int ch,int numsamples,BufIO &buf);
    int EncodeMonoFrame_Mapped(int ch,int numsamples,BufIO &buf);
    void Optimize(const FrameCoder::toptim_cfg &ocfg,SacProfile &profile,const std::vector<int>&params_to_optimize);
    double GetCost(const CostFunction *func,const tch_samples &samples,std::size_t samples_to_optimize) const;
    void PredictFrame(const SacProfile &profile,tch_samples &error,int from,int numsamples,bool optimize);
    void UnpredictFrame(const SacProfile &profile,int numsamples);
    double CalcRemapError(int ch, int numsamples);
    void EncodeMonoFrame(int ch,int numsamples);
    void DecodeMonoFrame(int ch,int numsamples);
    int numchannels_,framesize_,numsamples_;
    int profile_size_bytes_;
    SacProfile base_profile;
    coder_ctx opt;
};

class Codec {
  struct tsub_frame {
    int state=-1;
    int start=0;
    int length=0;
  };
  public:
    Codec(){};
    Codec(FrameCoder::coder_ctx &opt):opt_(opt) {};
    void EncodeFile(Wav &myWav,Sac &mySac);
    //void EncodeFile(Wav &myWav,Sac &mySac,int profile,int optimize,int sparse_pcm);
    void DecodeFile(Sac &mySac,Wav &myWav);
    void ScanFrames(Sac &mySac);
  private:
    std::vector<Codec::tsub_frame> Analyse(const std::vector <std::vector<int32_t>>&samples,int blocksamples,int min_frame_length,int samples_read);
    void PushState(std::vector<Codec::tsub_frame> &sub_frames,Codec::tsub_frame &curframe,int min_frame_length,int block_state,int samples_block);
    std::pair<double,double> AnalyseSparse(std::span<const int32_t> buf);
    void PrintProgress(int samplesprocessed,int totalsamples);
    FrameCoder::coder_ctx opt_;
    //int framesize;
};

#endif
