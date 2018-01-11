//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef WBWRITER_H
#define WBWRITER_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif


#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>


#include <basics/utilities/UbException.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>
#include <basics/utilities/UbPointerWrapper.h>
#include <basics/utilities/UbAutoRun.hpp>
#include <basics/objects/ObFactory.h>

class WbWriter
{
public:
   OBCREATOR_EXT(WbWriter)

   //////////////////////////////////////////////////////////////////////////
   virtual ~WbWriter() 
   {

   }

   //////////////////////////////////////////////////////////////////////////
   //rein virtuelle Methoden
   virtual std::string getFileExtension() = 0;

   //////////////////////////////////////////////////////////////////////////
   //nodes
   virtual std::string writeNodes(const std::string& filename,std::vector< UbTupleFloat3 >& nodes) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeNodesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeNodesWithNodeDataDouble(const std::string& filename,std::vector< UbTupleDouble3 >& nodes, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }

   //////////////////////////////////////////////////////////////////////////
   //lines
   //     0 ---- 1
   //nodenumbering must start with 0!
   virtual std::string writeLines(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeLinesWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }

   //////////////////////////////////////////////////////////////////////////
   //triangles
   //cell numbering:
   //                     2
   //                      
   //                  0 === 1
   //nodenumbering must start with 0!
   virtual std::string writeTriangles(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeTrianglesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }

   //////////////////////////////////////////////////////////////////////////
   //quads
   //cell numbering:
   //                  3---2
   //                  |   |
   //                  0---1
   //nodenumbering must start with 0!
   virtual std::string writeQuads(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeQuadsWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeQuadsWithCellData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >& datanames, std::vector< std::vector< double > >& celldata){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeQuadsWithNodeAndCellData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, 
                                                     std::vector< std::string >& nodedatanames, std::vector< std::vector< double > >& nodedata, std::vector< std::string >& celldatanames,
                                                     std::vector< std::vector< double > >& celldata                                                                                       ){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }

   //////////////////////////////////////////////////////////////////////////
   //octs
   //     7 ---- 6
   //    /|     /|
   //   4 +--- 5 |
   //   | |    | |
   //   | 3 ---+ 2
   //   |/     |/
   //   0 ---- 1
   virtual std::string writeOcts(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt8 >& cells){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeOctsWithCellData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& celldata){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }
   virtual std::string writeOctsWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata){ throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() );  }

private:

};


#ifdef CAB_RCF
//serialize von singletons muss hier etwas anders erfolgen ;-)
template<class Archive>
inline bool serializeWbWriter(Archive &ar, WbWriter*& writer)
{
   std::string writerID;

   if( ArchiveTools::isReading(ar) )
   {                                                                  
      ar & writerID;
      if(writerID!="no_WbWriter") writer = ObFactory<WbWriter>::getInstance()->createObject(writerID);
      else                        writer = NULL;
   }                                                                  
   else /* if (ar.isWrite())) if(Archive::is_saving())*/                                      
   {                                                                   
      if(writer) writerID = writer->getClassObjectTypeID(); 
      else       writerID = "no_WbWriter";
      ar & writerID;
   } 
   return true;
}
//////////////////
template<class Archive, class STL_container>
inline bool serializeWbWriter(Archive &ar, STL_container& writers)
{
   int       nofCounter;
   std::string    writerID;
   WbWriter* dummy;

   if( ArchiveTools::isReading(ar) )
   {                                                                  
      ar & nofCounter;
      for(int i=0; i<nofCounter; i++)
      {
         serializeWbWriter(ar, dummy);
         writers.push_back(dummy);
      }
   }                                                                  
   else                                 
   {                                                                   
      nofCounter = (int)writers.size();
      ar & nofCounter;
      typename STL_container::iterator pos;
      for(pos=writers.begin(); pos!=writers.end(); ++pos)
         serializeWbWriter(ar, *pos);
   }                                                                   

   return true;
}
//////////////////////////////////////////////////////////////////////////
// Spezialisierung des UbPointerWrappers fuer WbWriter... 
// da man bei singletons keine serializemethode einbauen kann...
template< >
class UbPointerWrapper< WbWriter > 
{
public:
   UbPointerWrapper() : pointer(NULL) {}

   UbPointerWrapper(WbWriter* pointer) : pointer(pointer) {}

   WbWriter* get() { return pointer; }

   template<class Archive>
   void serialize(Archive& ar, const unsigned int version) 
   {
      serializeWbWriter(ar, pointer);
   }

private:
   WbWriter* pointer;
};


#endif //CAB_RCF


#endif //WBWRITER_H
