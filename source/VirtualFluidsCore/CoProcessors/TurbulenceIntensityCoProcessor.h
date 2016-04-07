#ifndef TurbulenceIntensityCoProcessor_H
#define TurbulenceIntensityCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class TurbulenceIntensityCoProcessor;
typedef boost::shared_ptr<TurbulenceIntensityCoProcessor> TurbulenceIntensityCoProcessorPtr;

class TurbulenceIntensityCoProcessor : public CoProcessor
{
public:
   TurbulenceIntensityCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer, 
                          UbSchedulerPtr s, CommunicatorPtr comm);
   void process(double step);
protected:
   void collectData(double step);
   void addData(const Block3DPtr block);
   void clearData();
   void calculateAverageValues(double timeStep);
private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::vector<std::vector<Block3DPtr> > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   std::string path;
   WbWriter* writer;
   CommunicatorPtr comm;
   enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxxyyzz = 3};
};
#endif
