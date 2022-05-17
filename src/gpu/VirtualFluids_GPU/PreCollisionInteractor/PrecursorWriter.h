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
    uint nPoints;
    uint *indicesH, *indicesD;
    real *vxH, *vxD;
    real *vyH, *vyD;
    real *vzH, *vzD;
    UbTupleInt4 extent;
    UbTupleFloat3 origin;
    UbTupleFloat3 spacing;
    int* indicesOnPlane;
    uint ny, nz;
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
        uint _tWrite
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
    tWrite(_tWrite)
    {
        nFilesWritten = 0;
    };
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
private:
    WbWriterVtkXmlImageBinary* getWriter(){ return WbWriterVtkXmlImageBinary::getInstance(); };
    void write(Parameter* para, int level);
    void writeGridFile(Parameter* para, int level, uint part, int numberOfTimestepsPerPart);
    void writeParallelFile(Parameter* para, int level, uint parts);
    std::string makeGridFileName(int level, int id, uint part);
    std::string makeParallelFileName(int level, int id);

private:
    std::vector<SPtr<PrecursorStruct>> precursorStructs;
    std::vector<std::vector<std::vector<real>>> vx, vy, vz; // level, time, array
    std::string fileName, outputPath;
    std::vector<std::string> nodedatanames = {"vx", "vy", "vz"};
    std::vector<std::string> celldatanames, pieceSources;
    uint tStartOut, tSave, tWrite, nFilesWritten;
    real xPos, yMin, yMax, zMin, zMax;
    std::vector<UbTupleInt6> pieceExtents;
};

#endif //PRECURSORPROBE_H_