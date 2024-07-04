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

int Sac::WriteSACHeader(Wav &myWav)
{
  Chunks &myChunks=myWav.GetChunks();
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
  buf[16] = mcfg.profile;
  buf[17] = mcfg.max_framelen;

  // write wav meta data
  const uint32_t metadatasize=myChunks.GetMetaDataSize();
  BitUtils::put32LH(buf+18,metadatasize);
  file.write((char*)buf,22);
  if (myChunks.PackMetaData(metadata)!=metadatasize) std::cerr << "  warning: metadatasize mismatch\n";
  WriteData(metadata,metadatasize);
  return 0;
}

int Sac::UnpackMetaData(Wav &myWav)
{
  size_t unpackedbytes=myWav.GetChunks().UnpackMetaData(metadata);
  if (mcfg.metadatasize!=unpackedbytes) {std::cerr << "  warning: unpackmetadata mismatch\n";return 1;}
  return 0;
}

int Sac::ReadSACHeader()
{
  uint8_t buf[32];
  file.read((char*)buf,22);
  if (buf[0]=='S' && buf[1]=='A' && buf[2]=='C' && buf[3]=='2') {
    numchannels=BitUtils::get16LH(buf+4);
    samplerate=BitUtils::get32LH(buf+6);
    bitspersample=BitUtils::get16LH(buf+10);
    numsamples=BitUtils::get32LH(buf+12);
    mcfg.profile=buf[16];
    mcfg.max_framelen=buf[17];
    mcfg.metadatasize=BitUtils::get32LH(buf+18);
    ReadData(metadata,mcfg.metadatasize);
    mcfg.max_framesize=samplerate*static_cast<uint32_t>(mcfg.max_framelen);
    return 0;
  } else return 1;
}
