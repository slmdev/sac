#ifndef SAC_H
#define SAC_H

#include "file.h"
#include "wav.h"
#include "../common/utils.h"

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
    int WriteHeader(Wav &myWav);
    int ReadHeader();
    int UnpackMetaData(Wav &myWav);
    std::vector <uint8_t>metadata;
  private:
     uint32_t metadatasize,profile;
};


#endif // SAC_H
