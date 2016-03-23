#ifndef TurbulenceIntensityPostprocessor_H
#define TurbulenceIntensityPostprocessor_H

#include "Postprocessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class TurbulenceIntensityPostprocessor;
typedef boost::shared_ptr<TurbulenceIntensityPostprocessor> TurbulenceIntensityPostprocessorPtr;

class TurbulenceIntensityPostprocessor : public Postprocessor
{
public:
   TurbulenceIntensityPostprocessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer, 
                          UbSchedulerPtr s, CommunicatorPtr comm);
   void update(double step);
protected:
   void collectPostprocessData(double step);
   void addPostprocessData(const Block3DPtr block);
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
