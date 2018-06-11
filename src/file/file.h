#ifndef FILE_H
#define FILE_H

#include "../global.h"

class AudioFile {
  public:
    AudioFile():samplerate(0),bitspersample(0),numchannels(0),numsamples(0),kbps(0){};
    AudioFile(const AudioFile &file)
    :samplerate(file.getSampleRate()),bitspersample(file.getBitsPerSample()),
    numchannels(file.getNumChannels()),numsamples(file.getNumSamples()),kbps(0){};

    int OpenRead(const std::string &fname);
    int OpenWrite(const std::string &fname);
    std::streampos getFileSize() const {return filesize;};
    int getNumChannels()const {return numchannels;};
    int getSampleRate()const {return samplerate;};
    int getBitsPerSample()const {return bitspersample;};
    int getKBPS()const {return kbps;};
    int getNumSamples()const {return numsamples;};
    std::streampos readFileSize();
    void Close() {if (file.is_open()) file.close();};
    void ReadData(std::vector <uint8_t>&data,size_t len);
    void WriteData(const std::vector <uint8_t>&data,size_t len);
    std::fstream file;
  protected:
    std::streampos filesize;
    int samplerate,bitspersample,numchannels,numsamples,kbps;
};
#endif // FILE_H
