#ifndef Y_2D_Slice_To_RESULTS_H
#define Y_2D_Slice_To_RESULTS_H

#include "../ToVectorWriter.h"

#include <vector>
#include <memory>

class Parameter;

class Y2dSliceToResults : public ToVectorWriter
{
public:
	Y2dSliceToResults(std::shared_ptr<SimulationResults> simResults, unsigned int ySliceForCalculation, unsigned int startTimeY2dSliceToVector, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, std::shared_ptr<FileWriter> fileWriter, unsigned int startTimeDataWriter);

private:
	void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level);
	
	std::shared_ptr<SimulationResults> simResults;
	int CoordPara3DTo1D(int x, int y, int z);
	int CoordResults2DTo1D(int x, int z);
};
#endif