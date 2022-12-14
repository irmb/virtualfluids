#ifndef TRANSIENTBCSETTER_H_
#define TRANSIENTBCSETTER_H_

#include "Core/DataTypes.h"
#include <Core/StringUtilities/StringUtil.h>
#include "PointerDefinitions.h"

#include <string>
#include <vector>
#include <math.h>
#include <sstream>
#include <future>
class Grid;
namespace gg
{
    class BoundaryCondition;
}


enum class FileType
{
    VTK
};

struct Quantity
{
    std::string name;
    int offset;
    std::vector<double> values;
};

class VTKFile
{
public: 
    VTKFile(std::string _fileName): 
    fileName(_fileName)
    {
        readHeader();
        this->loaded = false;
        // printFileInfo();
    };

    void getData(real* data, uint numberOfNodes, const std::vector<uint>& readIndices, const std::vector<uint>& writeIndices, uint offsetRead, uint offsetWrite);
    bool markNANs(std::vector<uint> readIndices);
    bool inBoundingBox(real posX, real posY, real posZ){return  inXBounds(posX) && inYBounds(posY) && inZBounds(posZ); };
    bool inXBounds(real posX){ return posX<=maxX && posX>=minX; };
    bool inYBounds(real posY){ return posY<=maxY && posY>=minY; };
    bool inZBounds(real posZ){ return posZ<=maxZ && posZ>=minZ; };
    int findNeighborWSB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)  , getIdxSY(posY)  , getIdxBZ(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborWST(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)  , getIdxSY(posY)  , getIdxBZ(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborWNB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)  , getIdxSY(posY)+1, getIdxBZ(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborWNT(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)  , getIdxSY(posY)+1, getIdxBZ(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborESB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)+1, getIdxSY(posY)  , getIdxBZ(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborEST(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)+1, getIdxSY(posY)  , getIdxBZ(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborENB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)+1, getIdxSY(posY)+1, getIdxBZ(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborENT(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxWX(posX)+1, getIdxSY(posY)+1, getIdxBZ(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int getIdxX(int linearIdx){ return linearIdx%nx;};
    int getIdxY(int linearIdx){ return (linearIdx/nx)%ny;};
    int getIdxZ(int linearIdx){ return linearIdx/(nx*ny); };
    real getX(int linearIdx){ return getIdxX(linearIdx)*deltaX+minX; };
    real getY(int linearIdx){ return getIdxY(linearIdx)*deltaY+minY; };
    real getZ(int linearIdx){ return getIdxZ(linearIdx)*deltaZ+minZ; };
    int getIdxWX(real posX){ return (posX-minX)/deltaX; };
    int getIdxSY(real posY){ return (posY-minY)/deltaY; };
    int getIdxBZ(real posZ){ return (posZ-minZ)/deltaZ; };
    int getClosestIdxX(real posX){ int x = round((posX-minX)/deltaX); return x>nx ? nx : (x<0 ? 0 : x);};
    int getClosestIdxY(real posY){ int y = round((posY-minY)/deltaY); return y>ny ? ny : (y<0 ? 0 : y);};
    int getClosestIdxZ(real posZ){ int z = round((posZ-minZ)/deltaZ); return z>nz ? nz : (z<0 ? 0 : z);};
    int getLinearIndex(int idxX, int idxY, int idxZ){ return idxX + nx*(idxY+ny*idxZ); };
    int getNumberOfPointsInXYPlane(){ return nx*ny; }
    int getNumberOfPointsInYZPlane(){ return ny*nz; }
    int getNumberOfPointsInXZPlane(){ return nx*nz; }
    int getNumberOfPoints(){ return nx*ny*nz; }
    size_t getNumberOfQuantities(){ return quantities.size(); }
    void loadFile();
    void unloadFile();
    bool isLoaded(){return loaded;};


private:
    void readHeader();
    void printFileInfo();

public:

private:
    std::string fileName;
    real minX, maxX, minY, maxY, minZ, maxZ;
    real deltaX, deltaY, deltaZ;
    int nx, ny, nz;
    std::vector<Quantity> quantities;
    bool loaded;
};

class FileCollection
{
public:
    FileCollection(std::string _prefix): 
    prefix(_prefix){};

    virtual ~FileCollection() = default;

    virtual size_t getNumberOfQuantities()=0;

    virtual FileType getFileType()=0;

protected:
    std::string prefix;
};


class VTKFileCollection : public FileCollection
{
public:
    VTKFileCollection(std::string _prefix): 
    FileCollection(_prefix)
    {
        findFiles();
    };

    FileType getFileType(){ return FileType::VTK; };
    size_t getNumberOfQuantities(){ return files[0][0][0].getNumberOfQuantities(); }
    

private:
    void findFiles();
    std::string makeFileName(int level, int id, int part)
    { 
        return prefix + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" +    StringUtil::toString<int>(id)
                    + "_File_" +  StringUtil::toString<int>(part) 
                    + ".bin." + suffix;
    };


public:
    static const inline std::string suffix = "vti";
    std::vector<std::vector<std::vector<VTKFile>>> files;
};


class TransientBCInputFileReader
{
public:
    TransientBCInputFileReader()
    { 
        this->nPoints = 0; 
        this->nPointsRead = 0;
        this->writingOffset = 0;        
    };
    virtual ~TransientBCInputFileReader() = default;

    virtual void getNextData(real* data, uint numberOfNodes, real time)=0;
    virtual void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ)=0;
    uint getNPoints(){return nPoints; };
    uint getNPointsRead(){return nPointsRead; };
    size_t getNumberOfQuantities(){ return nQuantities; };
    void setWritingOffset(uint offset){ this->writingOffset = offset; }
    void getNeighbors(uint* neighborNT, uint* neighborNB, uint* neighborST, uint* neighborSN);
    void getWeights(real* _weightsNT, real* _weightsNB, real* _weightsST, real* _weightsSB);

public:
    std::vector<uint> planeNeighborNT,  planeNeighborNB, planeNeighborST, planeNeighborSB;
    std::vector<real> weightsNT, weightsNB, weightsST,  weightsSB;

protected:
    uint nPoints, nPointsRead, writingOffset;
    uint nReads=0;
    size_t nQuantities=0;
};


class VTKReader : public TransientBCInputFileReader
{
public:
    VTKReader(SPtr<VTKFileCollection> _fileCollection, uint _readLevel):
    fileCollection(_fileCollection), 
    readLevel(_readLevel)
    {
        this->nQuantities = fileCollection->getNumberOfQuantities();
        read = std::async([](){});
    };
    void getNextData(real* data, uint numberOfNodes, real time) override;
    void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ) override;
private:  
    uint getWriteIndex(int level, int id, int linearIdx);
    void initializeIndexVectors();

private:
    std::vector<std::vector<std::vector<uint>>> readIndices, writeIndices;
    std::vector<std::vector<size_t>> nFile;
    SPtr<VTKFileCollection> fileCollection;
    uint readLevel;
    std::future<void> read;
};


SPtr<FileCollection> createFileCollection(std::string prefix, FileType type);
SPtr<TransientBCInputFileReader> createReaderForCollection(SPtr<FileCollection> fileCollection, uint readLevel);

#endif //TRANSIENTBCSETTER_H_