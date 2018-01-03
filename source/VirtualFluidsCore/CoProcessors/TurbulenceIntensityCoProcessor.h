#ifndef TurbulenceIntensityCoProcessor_H
#define TurbulenceIntensityCoProcessor_H

#include <memory>
#include <string>
#include <string>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

class TurbulenceIntensityCoProcessor;
typedef std::shared_ptr<TurbulenceIntensityCoProcessor> TurbulenceIntensityCoProcessorPtr;

class TurbulenceIntensityCoProcessor : public CoProcessor
{
public:
   TurbulenceIntensityCoProcessor(std::shared_ptr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm);
   void process(double step);
protected:
   void collectData(double step);
   void addData(const std::shared_ptr<Block3D> block);
   void clearData();
   void calculateAverageValues(double timeStep);
private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::vector<std::vector<std::shared_ptr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   std::string path;
   WbWriter* writer;
   std::shared_ptr<Communicator> comm;
   enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxxyyzz = 3};
};
#endif
