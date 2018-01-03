#ifndef D3Q27MACROSCOPICQUANTITIESCoProcessor_H
#define D3Q27MACROSCOPICQUANTITIESCoProcessor_H

#include <memory>
#include <string>
#include <vector>

#include "CoProcessor.h"

#include "LBMSystem.h"

class Communicator;
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class WbWriter;
class Block3D;

class WriteMacroscopicQuantitiesCoProcessor;
typedef std::shared_ptr<WriteMacroscopicQuantitiesCoProcessor> MacroscopicQuantitiesCoProcessorPtr;

class WriteMacroscopicQuantitiesCoProcessor : public CoProcessor 
{
public:
   WriteMacroscopicQuantitiesCoProcessor();
   WriteMacroscopicQuantitiesCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s,
                                           const std::string& path, WbWriter* const writer, 
                                           std::shared_ptr<LBMUnitConverter> conv, std::shared_ptr<Communicator> comm);
   ~WriteMacroscopicQuantitiesCoProcessor(){}

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataMQ(std::shared_ptr<Block3D> block);
   void addDataGeo(std::shared_ptr<Block3D> block);
   void clearData();

private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::string path;
   WbWriter* writer;
   std::shared_ptr<LBMUnitConverter> conv;
   bool bcInformation;
   std::vector<std::vector<std::shared_ptr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   std::shared_ptr<Communicator> comm;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;

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
