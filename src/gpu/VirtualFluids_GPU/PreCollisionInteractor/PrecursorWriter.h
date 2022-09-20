#ifndef PRECURSORPROBE_H_
#define PRECURSORPROBE_H_

#include "PreCollisionInteractor.h"
#include "WbWriterVtkXmlImageBinary.h"
#include "LBM/LB.h"
#include <string>
#include <vector>
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

struct PrecursorStruct
{
    uint nPoints, nPointsInPlane, timestepsPerFile, filesWritten, timestepsBuffered;
    uint *indicesH, *indicesD;
    real *vxH, *vxD;
    real *vyH, *vyD;
    real *vzH, *vzD;
    DistributionSubset9 distH, distD;
    UbTupleInt4 extent;
    UbTupleFloat2 origin;
    UbTupleFloat3 spacing;
    int* indicesOnPlane;
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
        if(      _outputVariable==OutputVariable::Velocities   ){   nodedatanames = {"vx", "vy", "vz"};}
        else if( _outputVariable==OutputVariable::Distributions){   nodedatanames = {"f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8"};}
        else{ throw std::runtime_error("Invalid OutputVariable for PrecursorWriter"); }
    };
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    OutputVariable getOutputVariable(){ return this->outputVariable; }

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
    static std::string makeFileName(std::string fileName, int level, int id, uint part);
private:
    WbWriterVtkXmlImageBinary* getWriter(){ return WbWriterVtkXmlImageBinary::getInstance(); };
    void write(Parameter* para, int level);

private:
    std::vector<SPtr<PrecursorStruct>> precursorStructs;
    std::string fileName, outputPath;
    std::vector<std::string> nodedatanames;
    std::vector<std::string> celldatanames;
    uint tStartOut, tSave, maxtimestepsPerFile;
    real xPos, yMin, yMax, zMin, zMax;
    OutputVariable outputVariable;
};

#endif //PRECURSORPROBE_H_