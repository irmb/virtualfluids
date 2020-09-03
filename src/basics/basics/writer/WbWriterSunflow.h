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
private:
   WbWriterSunflow() : WbWriter() 
   {
      if(sizeof(unsigned char)!=1) throw UbException(UB_EXARGS,"error char  type mismatch");
      if(sizeof(int)          !=4) throw UbException(UB_EXARGS,"error int   type mismatch");
      if(sizeof(float)        !=4) throw UbException(UB_EXARGS,"error float type mismatch");
   }
   WbWriterSunflow( const WbWriterSunflow& );                  //no copy allowed 
   const WbWriterSunflow& operator=( const WbWriterSunflow& ); //no copy allowed

   static std::string  pvdEndTag;

public:
   std::string getFileExtension()  { return "ascii.sunflow"; }

   std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles);
};

#endif //WbWriterSunflow_H
