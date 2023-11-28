#ifndef TRANSIENTBCSETTER_H_
#define TRANSIENTBCSETTER_H_

#include <cmath>
#include <future>
#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <StringUtilities/StringUtil.h>


class Grid;
namespace gg
{
    class BoundaryCondition;
}

enum class TransientBCFileType
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
    explicit VTKFile(std::string fileName): 
    fileName(fileName)
    {
        readHeader();
        this->loaded = false;
        // printFileInfo();
    };

    void getData(real* data, uint numberOfNodes, const std::vector<uint>& readIndices, const std::vector<uint>& writeIndices, uint offsetRead, uint offsetWrite);
    bool markNANs(const std::vector<uint>& readIndices) const;
    bool inBoundingBox(real posX, real posY, real posZ) const{return  inXBounds(posX) && inYBounds(posY) && inZBounds(posZ); };
    bool inXBounds(real posX) const{ return posX<=maxX && posX>=minX; };
    bool inYBounds(real posY) const{ return posY<=maxY && posY>=minY; };
    bool inZBounds(real posZ) const{ return posZ<=maxZ && posZ>=minZ; };
    int findNeighborMMM(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)  , getIdx0M0(posY)  , getIdx00M(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborMMP(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)  , getIdx0M0(posY)  , getIdx00M(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborMPM(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)  , getIdx0M0(posY)+1, getIdx00M(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborMPP(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)  , getIdx0M0(posY)+1, getIdx00M(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborPMM(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)+1, getIdx0M0(posY)  , getIdx00M(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborPMP(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)+1, getIdx0M0(posY)  , getIdx00M(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborPPM(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)+1, getIdx0M0(posY)+1, getIdx00M(posZ)  ); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int findNeighborPPP(real posX, real posY, real posZ) const{ const int idx = getLinearIndex(getIdxM00(posX)+1, getIdx0M0(posY)+1, getIdx00M(posZ)+1); return (idx>=0) && (idx<nx*ny*nz) ? idx : -1; };
    int getIdxX(int linearIdx) const{ return linearIdx%nx;};
    int getIdxY(int linearIdx) const{ return (linearIdx/nx)%ny;};
    int getIdxZ(int linearIdx) const{ return linearIdx/(nx*ny); };
    real getX(int linearIdx) const{ return getIdxX(linearIdx)*deltaX+minX; };
    real getY(int linearIdx) const{ return getIdxY(linearIdx)*deltaY+minY; };
    real getZ(int linearIdx) const{ return getIdxZ(linearIdx)*deltaZ+minZ; };
    int getIdxM00(real posX) const{ return (posX-minX)/deltaX; };
    int getIdx0M0(real posY) const{ return (posY-minY)/deltaY; };
    int getIdx00M(real posZ) const{ return (posZ-minZ)/deltaZ; };
    int getClosestIdxX(real posX) const{ const int x = round((posX-minX)/deltaX); return x>nx ? nx : (x<0 ? 0 : x);};
    int getClosestIdxY(real posY) const{ const int y = round((posY-minY)/deltaY); return y>ny ? ny : (y<0 ? 0 : y);};
    int getClosestIdxZ(real posZ) const{ const int z = round((posZ-minZ)/deltaZ); return z>nz ? nz : (z<0 ? 0 : z);};
    int getLinearIndex(int idxX, int idxY, int idxZ) const{ return idxX + nx*(idxY+ny*idxZ); };
    int getNumberOfPointsInXYPlane() const{ return nx*ny; }
    int getNumberOfPointsInYZPlane() const{ return ny*nz; }
    int getNumberOfPointsInXZPlane() const{ return nx*nz; }
    int getNumberOfPoints() const{ return nx*ny*nz; }
    size_t getNumberOfQuantities(){ return quantities.size(); }
    void loadFile();
    void unloadFile();
    bool isLoaded() const{return loaded;};


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
    FileCollection(std::string prefix): 
    prefix(prefix){};

    virtual ~FileCollection() = default;

    virtual size_t getNumberOfQuantities() = 0;

    virtual TransientBCFileType getFileType() = 0;

protected:
    std::string prefix;
};


class VTKFileCollection : public FileCollection
{
public:
    VTKFileCollection(std::string prefix): 
    FileCollection(prefix)
    {
        findFiles();
    };

    TransientBCFileType getFileType() override{ return TransientBCFileType::VTK; };
    size_t getNumberOfQuantities() override{ return files[0][0][0].getNumberOfQuantities(); }
    

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
    uint getNPoints() const{return nPoints; };
    uint getNPointsRead() const{return nPointsRead; };
    size_t getNumberOfQuantities() const{ return nQuantities; };
    void setWritingOffset(uint offset){ this->writingOffset = offset; }
    void getNeighbors(uint* neighbor0PP, uint* neighbor0PM, uint* neighbor0MP, uint* neighbor0MM);
    void getWeights(real* _weights0PP, real* _weights0PM, real* _weights0MP, real* _weights0MM);

public:
    std::vector<uint> planeNeighbor0PP,  planeNeighbor0PM, planeNeighbor0MP, planeNeighbor0MM;
    std::vector<real> weights0PP, weights0PM, weights0MP,  weights0MM;

protected:
    uint nPoints, nPointsRead, writingOffset;
    uint nReads=0;
    size_t nQuantities=0;
};


class VTKReader : public TransientBCInputFileReader
{
public:
    VTKReader(SPtr<VTKFileCollection> fileCollection, uint readLevel):
    fileCollection(fileCollection), 
    readLevel(readLevel)
    {
        this->nQuantities = fileCollection->getNumberOfQuantities();
        read = std::async([](){});
    };
    void getNextData(real* data, uint numberOfNodes, real time) override;
    void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ) override;
private:  
    uint getWriteIndex(int level, int id, int linearIndex);
    void initializeIndexVectors();

private:
    std::vector<std::vector<std::vector<uint>>> readIndices, writeIndices;
    std::vector<std::vector<size_t>> nFile;
    SPtr<VTKFileCollection> fileCollection;
    uint readLevel;
    std::future<void> read;
};


SPtr<FileCollection> createFileCollection(std::string prefix, TransientBCFileType type);
SPtr<TransientBCInputFileReader> createReaderForCollection(SPtr<FileCollection> fileCollection, uint readLevel);

#endif //TRANSIENTBCSETTER_H_
