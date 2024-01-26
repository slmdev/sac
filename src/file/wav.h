#ifndef WAV_H
#define WAV_H

#include "file.h"
#include "../common/md5.h"

class Chunks {
  public:
    struct tChunk {
      uint32_t id,csize;
      std::vector <uint8_t>data;
    };
    Chunks():metadatasize(0){};
    void Append(uint32_t chunkid,uint32_t chunksize,const uint8_t *data,uint32_t len);
    size_t GetNumChunks() const {return wavchunks.size();};
    uint32_t GetChunkID(int chunk) const {return wavchunks[chunk].id;};
    uint32_t GetChunkSize(int chunk) const {return wavchunks[chunk].csize;};
    size_t GetChunkDataSize(int chunk) const {return wavchunks[chunk].data.size();};
    uint32_t GetMetaDataSize() const {return metadatasize;};
    size_t PackMetaData(std::vector <uint8_t>&data);
    size_t UnpackMetaData(const std::vector <uint8_t>&data);
    std::vector <tChunk> wavchunks;
  private:
    uint32_t metadatasize;
};

class Wav : public AudioFile {
  public:
    Wav(bool verbose=false);
    Wav(AudioFile &file,bool verbose=false);
    int ReadHeader();
    int WriteHeader();
    void InitFileBuf(int maxframesize);
    int ReadSamples(std::vector <std::vector <int32_t>>&data,int samplestoread);
    int WriteSamples(std::vector <std::vector <int32_t>>&data,int samplestowrite);
    Chunks &GetChunks(){return myChunks;};
    MD5::MD5Context md5ctx;
  private:
    Chunks myChunks;
    size_t chunkpos;
    std::vector <uint8_t>filebuffer;
    std::streampos datapos,endofdata;
    int byterate,blockalign,samplesleft;
    bool verbose;
};
#endif // WAV_H
