//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef WBWRITERVTKBINARY_H
#define WBWRITERVTKBINARY_H

#include <basics/writer/WbWriter.h>

#include <basics/memory/MbSharedPointerDefines.h>
class WbWriterVtkBinary;
typedef VFSharedPtr<WbWriterVtkBinary> WbWriterVtkBinaryPtr;


class WbWriterVtkBinary : public WbWriter
{
public:
   #ifndef SWIG 
   OBCREATOR_EXT( WbWriterVtkBinary )
   #endif

   static WbWriterVtkBinary* getInstance()
   {
      static WbWriterVtkBinary instance;
      return &instance;
   }
private:
   WbWriterVtkBinary() : WbWriter() {}                             
   WbWriterVtkBinary( const WbWriterVtkBinary& );                  //no copy allowed 
   const WbWriterVtkBinary& operator=( const WbWriterVtkBinary& ); //no copy allowed

public:
   std::string getFileExtension() { return ".bin.vtk"; }

   //////////////////////////////////////////////////////////////////////////
   //lines
   //     0 ---- 1
   //nodenumbering must start with 0!
   std::string writeLines(const std::string& filename, std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines);

   //////////////////////////////////////////////////////////////////////////
   //triangles
   //cell numbering:
   //                    2
   //                     
   //                  0---1
   //nodenumbering must start with 0!
   std::string writeTriangles(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells);
   std::string writeTrianglesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata);

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
   std::string writeOctsWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata);

};

#ifndef SWIG 
UB_AUTO_RUN_NAMED(ObFactory<WbWriter>::getInstance()->addObCreator(ObSingletonCreatorImpl<WbWriterVtkBinary ,WbWriter>::getInstance()), CAB_WbWriterVtkBinary);
#endif

#endif //WBWRITERVTKBINARY_H
