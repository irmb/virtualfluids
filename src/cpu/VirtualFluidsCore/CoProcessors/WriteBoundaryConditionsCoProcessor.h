#ifndef WriteBoundaryConditionsCoProcessor_H
#define WriteBoundaryConditionsCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "UbTuple.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;
class LBMUnitConverter;

class WriteBoundaryConditionsCoProcessor : public  CoProcessor
{
public:
   WriteBoundaryConditionsCoProcessor();
   WriteBoundaryConditionsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<Communicator> comm);
   ~WriteBoundaryConditionsCoProcessor() {}

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataGeo(SPtr<Block3D> block);
   void clearData();

private:
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleUInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::string path;
   WbWriter* writer;
   bool bcInformation;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   SPtr<Communicator> comm;
};
#endif
