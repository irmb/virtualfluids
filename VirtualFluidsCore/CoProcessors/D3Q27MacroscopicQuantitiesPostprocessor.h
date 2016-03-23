#ifndef D3Q27MACROSCOPICQUANTITIESPOSTPROCESSOR_H
#define D3Q27MACROSCOPICQUANTITIESPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class D3Q27MacroscopicQuantitiesPostprocessor;
typedef boost::shared_ptr<D3Q27MacroscopicQuantitiesPostprocessor> D3Q27MacroscopicQuantitiesPostprocessorPtr;

class D3Q27MacroscopicQuantitiesPostprocessor : public  Postprocessor 
{
public:
   D3Q27MacroscopicQuantitiesPostprocessor();
   D3Q27MacroscopicQuantitiesPostprocessor(Grid3DPtr grid, UbSchedulerPtr s, 
                                           const std::string& path, WbWriter* const writer, 
                                           LBMUnitConverterPtr conv,  
                                           bool nodesInformation = false);
   void update(double step);
protected:
   void collectPostprocessData(double step);
   void addPostprocessData1(Block3DPtr block);
   void addPostprocessData2(Block3DPtr block);
   void clearData();
private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::string path;
   WbWriter* writer;
   LBMUnitConverterPtr conv;
   bool bcInformation;
   std::vector<std::vector<Block3DPtr> > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<Postprocessor>(*this);
      ar & path;
      ar & conv;
      ar & bcInformation;
      ar & blockVector;
      ar & minInitLevel;
      ar & maxInitLevel;
      ar & gridRank;
      ar & writer;
   }
};
#endif
