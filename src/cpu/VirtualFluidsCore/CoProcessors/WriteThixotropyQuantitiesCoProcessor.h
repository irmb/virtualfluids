#ifndef WriteThixotropyQuantitiesCoProcessor_H
#define WriteThixotropyQuantitiesCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "WbWriter.h"

class WriteThixotropyQuantitiesCoProcessor : public  CoProcessor
{
public:
	WriteThixotropyQuantitiesCoProcessor();
	WriteThixotropyQuantitiesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<LBMUnitConverter> conv, SPtr<Communicator> comm);
	~WriteThixotropyQuantitiesCoProcessor() = default;

   void process(double step) override;

protected:
   void collectData(double step);
   void addDataMQ(SPtr<Block3D> block);
   void clearData();

private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleUInt8> cells;
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
//	double ConcentrationSum;
};
#endif
