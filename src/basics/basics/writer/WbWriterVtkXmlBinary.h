#ifndef WBWRITERVTKXMLBINARY_H
#define WBWRITERVTKXMLBINARY_H

#include <string>

#include "basics_export.h"

#include <basics/writer/WbWriter.h>

class BASICS_EXPORT WbWriterVtkXmlBinary  : public WbWriter
{
public:
#ifndef SWIG
   OBCREATOR_EXT( WbWriterVtkXmlBinary )
#endif

   static WbWriterVtkXmlBinary* getInstance()
   {
      static WbWriterVtkXmlBinary instance;
      return &instance;
   }
private:
   WbWriterVtkXmlBinary() : WbWriter() 
   {
      if(sizeof(unsigned char)!=1) throw UbException(UB_EXARGS,"machine error char  type mismatch");
      if(sizeof(int)          !=4) throw UbException(UB_EXARGS,"machine error int   type mismatch");
      if(sizeof(float)        !=4) throw UbException(UB_EXARGS,"machine error float type mismatch");
   }

   WbWriterVtkXmlBinary( const WbWriterVtkXmlBinary& );                  //no copy allowed 
   const WbWriterVtkXmlBinary& operator=( const WbWriterVtkXmlBinary& ); //no copy allowed

   static std::string  pvdEndTag;
public:
   std::string getFileExtension() { return ".bin.vtu";   }

   std::string getHeaderTag();

   //write a metafile 
   std::string writeCollection(const std::string& filename, const std::vector<std::string>& filenames, const double& timestep, const bool& sepGroups);
   std::string addFilesToCollection(const std::string& filename, const std::vector<std::string>& filenames, const double& timestep, const bool& sepGroups);
   std::string writeParallelFile(const std::string& filename,std::vector<std::string>& pieceSources, std::vector<std::string>& pointDataNames, std::vector<std::string>& cellDataNames);

   //////////////////////////////////////////////////////////////////////////
   //nodes
   std::string writeNodes(const std::string& filename,std::vector< UbTupleFloat3 >& nodes);
   std::string writeNodesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata);

   //////////////////////////////////////////////////////////////////////////
   //lines
   //     0 ---- 1
   //nodenumbering must start with 0!
   std::string writeLines(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines);
   std::string writeLinesWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata);

   //////////////////////////////////////////////////////////////////////////
   //triangles
   //                    2
   //                     
   //                  0---1
   //nodenumbering must start with 0!
   std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles);
   std::string writeTrianglesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata);

   //////////////////////////////////////////////////////////////////////////
   //2D
   //cell numbering:
   //                  3---2
   //                  |   |
   //                  0---1
   //nodenumbering must start with 0!

   std::string writeQuads(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells);
   std::string writeQuadsWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata);
   std::string writeQuadsWithCellData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >& datanames, std::vector< std::vector< double > >& celldata);
   std::string writeQuadsWithNodeAndCellData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, 
                                             std::vector< std::string >& nodedatanames, std::vector< std::vector< double > >& nodedata, std::vector< std::string >& celldatanames,
                                             std::vector< std::vector< double > >& celldata                                                                    );
   
   //////////////////////////////////////////////////////////////////////////
   //octs
   //     7 ---- 6
   //    /|     /|
   //   4 +--- 5 |
   //   | |    | |
   //   | 3 ---+ 2
   //   |/     |/
   //   0 ---- 1
   std::string writeOcts(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt8 >& cells);
   std::string writeOctsWithCellData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& celldata);
   std::string writeOctsWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleUInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata);
   
private:

};

#ifndef SWIG
UB_AUTO_RUN_NAMED(ObFactory<WbWriter>::getInstance()->addObCreator(ObSingletonCreatorImpl<WbWriterVtkXmlBinary ,WbWriter>::getInstance()), CAB_WbWriterVtkXmlBinary);
#endif

#endif //WBWRITERVTKXMLBINARY_H
