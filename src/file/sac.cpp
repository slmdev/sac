#include "sac.h"
#include "../common/utils.h"
#include <iostream>

std::streamsize Sac::WriteMD5(uint8_t digest[16])
{
  file.write(reinterpret_cast<char*>(digest),16);
  return file.gcount();
}

std::streamsize Sac::ReadMD5(uint8_t digest[16])
{
  file.read(reinterpret_cast<char*>(digest), 16);
  return file.gcount();
}

int Sac::WriteHeader(Wav &myWav)
{
  Chunks &myChunks=myWav.GetChunks();
  const uint32_t metadatasize=myChunks.GetMetaDataSize();
  uint8_t buf[32];
  std::vector <uint8_t>metadata;
  buf[0]='S';
  buf[1]='A';
  buf[2]='C';
  buf[3]='2';
  BitUtils::put16LH(buf+4,numchannels);
  BitUtils::put32LH(buf+6,samplerate);
  BitUtils::put16LH(buf+10,bitspersample);
  BitUtils::put32LH(buf+12,numsamples);
  BitUtils::put32LH(buf+16,profile);
  BitUtils::put32LH(buf+20,metadatasize);
  file.write((char*)buf,24);
  if (myChunks.PackMetaData(metadata)!=metadatasize) std::cerr << "  warning: metadatasize mismatch\n";
  WriteData(metadata,metadatasize);
  return 0;
}

int Sac::UnpackMetaData(Wav &myWav)
{
  size_t unpackedbytes=myWav.GetChunks().UnpackMetaData(metadata);
  if (metadatasize!=unpackedbytes) {std::cerr << "  warning: unpackmetadata mismatch\n";return 1;}
  else return 0;
}

int Sac::ReadHeader()
{
  uint8_t buf[32];
  file.read((char*)buf,24);
  if (buf[0]=='S' && buf[1]=='A' && buf[2]=='C' && buf[3]=='2') {
    numchannels=BitUtils::get16LH(buf+4);
    samplerate=BitUtils::get32LH(buf+6);
    bitspersample=BitUtils::get16LH(buf+10);
    numsamples=BitUtils::get32LH(buf+12);
    profile=BitUtils::get32LH(buf+16);
    metadatasize=BitUtils::get32LH(buf+20);
    ReadData(metadata,metadatasize);

    return 0;
  } else return 1;
}
