#ifndef TurbulenceIntensityCoProcessor_H
#define TurbulenceIntensityCoProcessor_H

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

class TurbulenceIntensityCoProcessor : public CoProcessor
{
public:
   TurbulenceIntensityCoProcessor(SPtr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       SPtr<UbScheduler> s, SPtr<Communicator> comm);
   void process(double step) override;
protected:
   void collectData(double step);
   void addData(const SPtr<Block3D> block);
   void clearData();
   void calculateAverageValues(double timeStep);
private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleUInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   std::string path;
   WbWriter* writer;
   SPtr<Communicator> comm;
   enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxxyyzz = 3};
};
#endif
