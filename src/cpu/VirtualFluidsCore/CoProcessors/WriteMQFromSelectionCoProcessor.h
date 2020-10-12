#ifndef WriteMQFromSelectionCoProcessor_H
#define WriteMQFromSelectionCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"

#include "LBMSystem.h"
#include "UbTuple.h"


class Communicator;
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class WbWriter;
class Block3D;
class GbObject3D;

class WriteMQFromSelectionCoProcessor : public CoProcessor 
{
public:
   WriteMQFromSelectionCoProcessor();
   WriteMQFromSelectionCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                           SPtr<GbObject3D> gbObject,
                                           const std::string& path, WbWriter* const writer, 
                                           SPtr<LBMUnitConverter> conv, SPtr<Communicator> comm);
   ~WriteMQFromSelectionCoProcessor() override= default;

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataMQ(SPtr<Block3D> block);
   void clearData();

private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::string path;
   WbWriter* writer;
   SPtr<LBMUnitConverter> conv;
//   bool bcInformation;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   SPtr<Communicator> comm;
   SPtr<GbObject3D> gbObject;

   using CalcMacrosFct = void (*)(const LBMReal *const &, LBMReal &, LBMReal &, LBMReal &, LBMReal &);
   CalcMacrosFct calcMacros;
};

#endif
