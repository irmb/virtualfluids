#ifndef PRECURSORPROBE_H_
#define PRECURSORPROBE_H_

#include "PreCollisionInteractor.h"
#include "WbWriterVtkXmlImageBinary.h"

#include <string>
#include <vector>
#include "PointerDefinitions.h"


class Parameter;
class CudaMemoryManager;
class GridProvider;

struct PrecursorStruct
{
    uint nPoints, nPointsInPlane, timestepsPerFile, filesWritten, timestepsBuffered;
    uint *indicesH, *indicesD;
    real *vxH, *vxD;
    real *vyH, *vyD;
    real *vzH, *vzD;
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
    maxtimestepsPerFile(_maxTimestepsPerFile)
    {};
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
    static std::string makeFileName(std::string fileName, int level, int id, uint part);
private:
    WbWriterVtkXmlImageBinary* getWriter(){ return WbWriterVtkXmlImageBinary::getInstance(); };
    void write(Parameter* para, int level);

private:
    std::vector<SPtr<PrecursorStruct>> precursorStructs;
    std::string fileName, outputPath;
    std::vector<std::string> nodedatanames = {"vx", "vy", "vz"};
    std::vector<std::string> celldatanames;
    uint tStartOut, tSave, maxtimestepsPerFile;
    real xPos, yMin, yMax, zMin, zMax;
};

#endif //PRECURSORPROBE_H_