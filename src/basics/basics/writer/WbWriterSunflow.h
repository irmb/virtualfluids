#ifndef WbWriterSunflow_H
#define WbWriterSunflow_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterSunflow  : public WbWriter
{
public:
   static WbWriterSunflow* getInstance()
   {
      static WbWriterSunflow instance;
      return &instance;
   }

    WbWriterSunflow( const WbWriterSunflow& ) = delete;
    const WbWriterSunflow& operator=( const WbWriterSunflow& ) = delete;
private:
   WbWriterSunflow() : WbWriter() 
   {
      if(sizeof(unsigned char)!=1) throw UbException(UB_EXARGS,"error char  type mismatch");
      if(sizeof(int)          !=4) throw UbException(UB_EXARGS,"error int   type mismatch");
      if(sizeof(float)        !=4) throw UbException(UB_EXARGS,"error float type mismatch");
   }

   static std::string  pvdEndTag;

public:
   std::string getFileExtension() override  { return "ascii.sunflow"; }

   std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles) override;
};

#endif //WbWriterSunflow_H
