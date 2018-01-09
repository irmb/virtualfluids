#ifndef WriteBoundaryConditionsCoProcessor_H
#define WriteBoundaryConditionsCoProcessor_H

#include <memory>
#include <string>
#include <vector>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;
class LBMUnitConverter;

class WriteBoundaryConditionsCoProcessor;
typedef std::shared_ptr<WriteBoundaryConditionsCoProcessor> WriteBoundaryConditionsCoProcessorPtr;

class WriteBoundaryConditionsCoProcessor : public  CoProcessor
{
public:
   WriteBoundaryConditionsCoProcessor();
   WriteBoundaryConditionsCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s,
      const std::string& path, WbWriter* const writer,
       std::shared_ptr<LBMUnitConverter> conv, std::shared_ptr<Communicator> comm);
   ~WriteBoundaryConditionsCoProcessor() {}

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataGeo(std::shared_ptr<Block3D> block);
   void clearData();

private:
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

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<CoProcessor>(*this);
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
