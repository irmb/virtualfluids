#ifdef CAB_ZLIB
   #ifndef WBWRITERBOBJ_H
   #define WBWRITERBOBJ_H

   #include <string>
   #include <basics/writer/WbWriter.h>

   class WbWriterBOBJ  : public WbWriter
   {
   public:
      OBCREATOR_EXT( WbWriterBOBJ )

         static WbWriterBOBJ* getInstance()
      {
         static WbWriterBOBJ instance;
         return &instance;
      }
   private:
      WbWriterBOBJ() : WbWriter() 
      {
         if(sizeof(unsigned char)!=1) throw UbException(UB_EXARGS,"error char  type mismatch");
         if(sizeof(int)          !=4) throw UbException(UB_EXARGS,"error int   type mismatch");
         if(sizeof(float)        !=4) throw UbException(UB_EXARGS,"error float type mismatch");
      }
      WbWriterBOBJ( const WbWriterBOBJ& );                  //no copy allowed 
      const WbWriterBOBJ& operator=( const WbWriterBOBJ& ); //no copy allowed

      static std::string  pvdEndTag;

   public:
      std::string getFileExtension()  { return "BOBJ.gz"; }

      std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles);
   };

   UB_AUTO_RUN_NAMED(ObFactory<WbWriter>::getInstance()->addObCreator(ObSingletonCreatorImpl<WbWriterBOBJ ,WbWriter>::getInstance()), CAB_WbWriterVtkXmlASCII);

   #endif //WBWRITERBOBJ_H

#endif //CAB_ZLIB
