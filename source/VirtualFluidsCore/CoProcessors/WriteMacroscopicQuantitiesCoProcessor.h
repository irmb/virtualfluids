#ifndef WriteMacroscopicQuantitiesCoProcessor_H
#define WriteMacroscopicQuantitiesCoProcessor_H

#include <PointerDefinitions.h>
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

class WriteMacroscopicQuantitiesCoProcessor : public CoProcessor 
{
public:
   WriteMacroscopicQuantitiesCoProcessor();
   WriteMacroscopicQuantitiesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                           const std::string& path, WbWriter* const writer, 
                                           SPtr<LBMUnitConverter> conv, SPtr<Communicator> comm);
   ~WriteMacroscopicQuantitiesCoProcessor(){}

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataMQ(SPtr<Block3D> block);
   void clearData();

private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data; 
   std::string path;
   WbWriter* writer;
   SPtr<LBMUnitConverter> conv;
   bool bcInformation;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   SPtr<Communicator> comm;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;
};

#endif
