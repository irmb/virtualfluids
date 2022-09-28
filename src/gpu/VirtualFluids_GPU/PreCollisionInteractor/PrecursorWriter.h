#ifndef PRECURSORPROBE_H_
#define PRECURSORPROBE_H_

#include "PreCollisionInteractor.h"
#include "WbWriterVtkXmlImageBinary.h"
#include "LBM/LB.h"
#include <string>
#include <vector>
#include <future>
#include "PointerDefinitions.h"


class Parameter;
class CudaMemoryManager;
class GridProvider;

enum class OutputVariable {
   //! - Velocities
    Velocities,
    //! - Distributions
    Distributions    
};

static constexpr uint PrecP00 = 0;
static constexpr uint PrecPP0 = 1;
static constexpr uint PrecPM0 = 2;
static constexpr uint PrecP0P = 3;
static constexpr uint PrecP0M = 4;
static constexpr uint PrecPPP = 5;
static constexpr uint PrecPMP = 6;
static constexpr uint PrecPPM = 7;
static constexpr uint PrecPMM = 8;

struct PrecursorStruct
{
    uint nPoints, nPointsInPlane, timestepsPerFile, filesWritten, timestepsBuffered;
    uint *indicesH, *indicesD;
    real *dataH, *dataD;
    real *bufferH, *bufferD;
    uint nQuantities;
    UbTupleInt4 extent;
    UbTupleFloat2 origin;
    UbTupleFloat3 spacing;
    int* indicesOnPlane;
    cudaStream_t stream;
};

class PrecursorWriter : public PreCollisionInteractor
{
public:
    PrecursorWriter(
        const std::string _fileName,
        const std::string _outputPath,
        real _xPos,
        real _yMin, real _yMax,
        real _zMin, real _zMax,
        uint _tStartOut,
        uint _tSave,
        OutputVariable _outputVariable,
        uint _maxTimestepsPerFile=uint(1e4)
    ): 
    fileName(_fileName), 
    outputPath(_outputPath), 
    xPos(_xPos),
    yMin(_yMin),
    yMax(_yMax),
    zMin(_zMin),
    zMax(_zMax),
    tStartOut(_tStartOut), 
    tSave(_tSave),
    outputVariable(_outputVariable),
    maxtimestepsPerFile(_maxTimestepsPerFile)
    {
        nodedatanames = determineNodeDataNames();
        writeFuture = std::async([](){});
    };

    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    OutputVariable getOutputVariable(){ return this->outputVariable; }

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
    static std::string makeFileName(std::string fileName, int level, int id, uint part);
    
private:
    WbWriterVtkXmlImageBinary* getWriter(){ return WbWriterVtkXmlImageBinary::getInstance(); };
    void write(Parameter* para, int level, int timestepsBuffered);

    std::vector<std::string> determineNodeDataNames()
    {
        switch (outputVariable)
        {
        case OutputVariable::Velocities:
            return {"vx", "vy", "vz"};
            break;       
        case OutputVariable::Distributions:
            return {"fP00", "fPP0", "fPM0", "fP0P", "fP0M", "fPPP", "fPMP", "fPPM", "fPMM"};
            break;
        
        default:
            throw std::runtime_error("Invalid OutputVariable for PrecursorWriter");
            break;
        }
    }

private:
    std::vector<SPtr<PrecursorStruct>> precursorStructs;
    std::string fileName, outputPath;
    std::vector<std::string> nodedatanames;
    std::vector<std::string> celldatanames;
    uint tStartOut, tSave, maxtimestepsPerFile;
    real xPos, yMin, yMax, zMin, zMax;
    OutputVariable outputVariable;
    std::future<void> writeFuture;
};

#endif //PRECURSORPROBE_H_