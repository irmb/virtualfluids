#ifndef WBWRITERX3D_H
#define WBWRITERX3D_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterX3D  : public WbWriter
{
public:
   OBCREATOR_EXT( WbWriterX3D )

   static WbWriterX3D* getInstance()
   {
      static WbWriterX3D instance;
      return &instance;
   }
private:
   WbWriterX3D() : WbWriter() 
   {
      if(sizeof(unsigned char)!=1) throw UbException(UB_EXARGS,"error char  type mismatch");
      if(sizeof(int)          !=4) throw UbException(UB_EXARGS,"error int   type mismatch");
      if(sizeof(float)        !=4) throw UbException(UB_EXARGS,"error float type mismatch");
   }
   WbWriterX3D( const WbWriterX3D& );                  //no copy allowed 
   const WbWriterX3D& operator=( const WbWriterX3D& ); //no copy allowed

   static std::string  pvdEndTag;

public:
   std::string getFileExtension()  { return "ascii.X3D"; }

   std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles);
};

UB_AUTO_RUN_NAMED(ObFactory<WbWriter>::getInstance()->addObCreator(ObSingletonCreatorImpl<WbWriterX3D ,WbWriter>::getInstance()), CAB_WbWriterX3D);

#endif //WBWRITERX3D_H
