#include "file.h"

std::streampos AudioFile::readFileSize()
{
    std::streampos oldpos=file.tellg();
    file.seekg(0,std::ios_base::end);
    std::streampos fsize = file.tellg();
    file.seekg(oldpos);
    return fsize;
}

int AudioFile::OpenRead(const std::string &fname)
{
    file.open(fname,std::ios_base::in|std::ios_base::binary);
    if (file.is_open()) {filesize=readFileSize();return 0;}
    else return 1;
}

int AudioFile::OpenWrite(const std::string &fname)
{
  file.open(fname,std::ios_base::out|std::ios_base::binary);
  if (file.is_open()) return 0;
  else return 1;
}

void AudioFile::ReadData(std::vector <uint8_t>&data,size_t len)
{
  if (data.size()<len) data.resize(len);
  file.read(reinterpret_cast<char*>(&data[0]),len);
}

void AudioFile::WriteData(const std::vector <uint8_t>&data,size_t len)
{
  file.write(reinterpret_cast<const char*>(&data[0]),len);
}
