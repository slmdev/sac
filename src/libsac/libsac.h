#ifndef CODEC_H
#define CODEC_H

#include "../file/wav.h"
#include "../file/sac.h"
#include "cost.h"
#include "profile.h"

class FrameCoder {
  public:
    struct coder_ctx {
      enum SearchCost {L1,Entropy,Golomb,Bitplane};
      enum SearchMethod {CMA,DDS};
      int optimize=0;
      int sparse_pcm=0;
      int optimize_maxnfunc=0;
      int verbose_level=0;
      int reset_profile=0;
      int zero_mean=0;
      int stereo_ms=0;
      int max_framelen=0;
      SearchMethod optimize_search=DDS;
      SearchCost optimize_cost=L1;
      double optimize_fraction=0;
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
    std::vector <std::vector<int32_t>>samples,err0,err1,error,s2u_error,s2u_error_map,pred;
    std::vector <BufIO> encoded,enc_temp1,enc_temp2;
    std::vector <SacProfile::FrameStats> framestats;

    static int WriteBlockHeader(std::fstream &file, const std::vector<SacProfile::FrameStats> &framestats, int ch);
    static int ReadBlockHeader(std::fstream &file, std::vector<SacProfile::FrameStats> &framestats, int ch);
  private:
    void EncodeProfile(const SacProfile &profile,std::vector <uint8_t>&buf);
    void DecodeProfile(SacProfile &profile,const std::vector <uint8_t>&buf);
    void AnalyseMonoChannel(int ch, int numsamples);
    void AnalyseShift(int ch,int numsamples);
    double AnalyseStereoChannel(int ch0, int ch1, int numsamples);
    void ApplyMs(int ch0, int ch1, int numsamples);
    //void InterChannel(int ch0,int ch1,int numsamples);
    int EncodeMonoFrame_Normal(int ch,int numsamples,BufIO &buf);
    int EncodeMonoFrame_Mapped(int ch,int numsamples,BufIO &buf);
    void Optimize(SacProfile &profile,const std::vector<int>&params_to_optimize);
    double GetCost(SacProfile &profile,CostFunction *func,int start_sample,int samples_to_optimize);

    void PredictFrame(const SacProfile &profile,int from,int numsamples,bool optimize=false);
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
  public:
    Codec(){};
    void EncodeFile(Wav &myWav,Sac &mySac,FrameCoder::coder_ctx &opt);
    //void EncodeFile(Wav &myWav,Sac &mySac,int profile,int optimize,int sparse_pcm);
    void DecodeFile(Sac &mySac,Wav &myWav);
    void ScanFrames(Sac &mySac);
  private:
    void PrintProgress(int samplesprocessed,int totalsamples);
    //int framesize;
};

#endif
