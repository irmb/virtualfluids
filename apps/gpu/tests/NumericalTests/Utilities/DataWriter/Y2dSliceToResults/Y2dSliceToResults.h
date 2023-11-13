#ifndef Y_2D_Slice_To_RESULTS_H
#define Y_2D_Slice_To_RESULTS_H

#include "../ToVectorWriter.h"

#include <vector>
#include <memory>

class Parameter;
class SimulationResults;
struct VectorWriterInformationStruct;

class Y2dSliceToResults : public ToVectorWriter
{
public:
    static std::shared_ptr<Y2dSliceToResults> getNewInstance(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults, unsigned int ySliceForCalculation);

    
private:
    Y2dSliceToResults();
    Y2dSliceToResults(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults, unsigned int ySliceForCalculation);
    
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level);
    int CoordPara3DTo1D(int x, int y, int z);
    int CoordResults2DTo1D(int x, int z);

    std::shared_ptr<SimulationResults> simResults;
    unsigned int ySliceForCalculation;
    unsigned int maxX, maxY, maxZ;
};
#endif