#ifndef D3Q27MACROSCOPICQUANTITIESCoProcessor_H
#define D3Q27MACROSCOPICQUANTITIESCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class WriteMacroscopicQuantitiesCoProcessor;
typedef boost::shared_ptr<WriteMacroscopicQuantitiesCoProcessor> MacroscopicQuantitiesCoProcessorPtr;

class WriteMacroscopicQuantitiesCoProcessor : public  CoProcessor 
{
public:
   WriteMacroscopicQuantitiesCoProcessor();
   WriteMacroscopicQuantitiesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, 
                                           const std::string& path, WbWriter* const writer, 
                                           LBMUnitConverterPtr conv,  
                                           CommunicatorPtr comm);
   ~WriteMacroscopicQuantitiesCoProcessor(){}
   void process(double step);
protected:
   void collectData(double step);
   void addDataMQ(Block3DPtr block);
   void addDataGeo(Block3DPtr block);
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
   CommunicatorPtr comm;

   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<CoProcessor>(*this);
   //   ar & path;
   //   ar & conv;
   //   ar & blockVector;
   //   ar & minInitLevel;
   //   ar & maxInitLevel;
   //   ar & gridRank;
   //   ar & writer;
   //}
};
#endif
