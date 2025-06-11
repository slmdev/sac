#ifndef SAC_H
#define SAC_H

#include "file.h"
#include "wav.h"

struct tFrameHeader {

};

class Sac : public AudioFile
{
  public:
    struct sac_cfg
    {
      uint8_t max_framelen=0;

      uint32_t max_framesize=0;
      uint32_t metadatasize=0;
    } mcfg;
    Sac():AudioFile(){};
    Sac(Wav &file)
    :AudioFile(file)
    {

    }
    void WriteFrameHeader(tFrameHeader &hdr);
    int WriteSACHeader(Wav &myWav);
    void WriteMD5(uint8_t digest[16]);
    void ReadMD5(uint8_t digest[16]);
    int ReadSACHeader();
    int UnpackMetaData(Wav &myWav);
    std::vector <uint8_t>metadata;
};


#endif // SAC_H
