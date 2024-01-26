#ifndef SAC_H
#define SAC_H

#include "file.h"
#include "wav.h"

struct tFrameHeader {

};

class Sac : public AudioFile
{
  public:
    Sac():AudioFile(),metadatasize(0),profile(0){};
    Sac(Wav &file)
    :AudioFile(file),metadatasize(0)
    {

    }
    void SetProfile(int val){profile=val;};
    int GetProfile(){return profile;};
    void WriteFrameHeader(tFrameHeader &hdr);
    int WriteHeader(Wav &myWav);
    std::streamsize WriteMD5(uint8_t digest[16]);
    std::streamsize ReadMD5(uint8_t digest[16]);
    int ReadHeader();
    int UnpackMetaData(Wav &myWav);
    std::vector <uint8_t>metadata;
  private:
     uint32_t metadatasize,profile;
};


#endif // SAC_H
