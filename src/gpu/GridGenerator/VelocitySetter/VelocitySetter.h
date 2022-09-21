#ifndef VELOCITY_SETTER_H_
#define VELOCITY_SETTER_H_

#include "Core/DataTypes.h"
#include <Core/StringUtilities/StringUtil.h>
#include "PointerDefinitions.h"

#include <string>
#include <vector>
#include <math.h>
#include <sstream>


class Grid;
namespace gg
{
    class BoundaryCondition;
}


enum class FileType
{
    VTK
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

    void getVelocities(real* vx, real* vy, real* vz, std::vector<uint> readIndeces, std::vector<uint> writeIndices, uint offsetRead, uint offsetWrite);
    void getDistributions(  real* f0, real* f1, real* f2, real* f3, real* f4, real* f5, real* f6, real* f7, real* f8, 
                            std::vector<uint> readIndeces, std::vector<uint> writeIndices, uint offsetRead, uint offsetWrite);
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
    int getClosestIdxX(real posX){ return round((posX-minX)/deltaX); };
    int getClosestIdxY(real posY){ return round((posY-minY)/deltaY); };
    int getClosestIdxZ(real posZ){ return round((posZ-minZ)/deltaZ); };
    int getLinearIndex(int idxX, int idxY, int idxZ){ return idxX + nx*(idxY+ny*idxZ); };
    int getNumberOfPointsInXYPlane(){ return nx*ny; }
    int getNumberOfPointsInYZPlane(){ return ny*nz; }
    int getNumberOfPointsInXZPlane(){ return nx*nz; }
    int getNumberOfPoints(){ return nx*ny*nz; }
    void loadFile();
    void unloadFile();


private:
    void readHeader();
    void printFileInfo();

public:

private:
    int offsetVx, offsetVy, offsetVz;
    std::string fileName;
    real minX, maxX, minY, maxY, minZ, maxZ;
    real deltaX, deltaY, deltaZ;
    int nx, ny, nz;
    std::vector<double> vxFile, vyFile, vzFile;
    std::vector<double> f0File, f1File, f2File, f3File, f4File, f5File, f6File, f7File, f8File; 
    bool loaded;
};


class VelocityFileCollection
{
public:
    VelocityFileCollection(std::string _prefix): 
    prefix(_prefix){};

    virtual ~VelocityFileCollection() = default;

    virtual FileType getFileType()=0;

protected:
    std::string prefix;
};


class VTKFileCollection : public VelocityFileCollection
{
public:
    VTKFileCollection(std::string _prefix): 
    VelocityFileCollection(_prefix)
    {
        findFiles();
    };

    FileType getFileType(){ return FileType::VTK; };
    

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


class VelocityReader
{
public:
    VelocityReader()
    { 
        this->nPoints = 0; 
        this->nPointsRead = 0;
        this->writingOffset = 0;
        
    };
    virtual ~VelocityReader() = default;

    virtual void getNextVelocities(real* vx, real* vy, real* vz, real t)=0;
    virtual void getNextDistributions(real* f0, real* f1, real* f2, real* f3, real* f4, real* f5, real* f6, real* f7, real* f8, real t)=0;
    virtual void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ)=0;
    uint getNPoints(){return nPoints; };
    uint getNPointsRead(){return nPointsRead; };
    void setWritingOffset(uint offset){ this->writingOffset = offset; }
    void getNeighbors(uint* neighborNT, uint* neighborNB, uint* neighborST, uint* neighborSN);
    void getWeights(real* _weightsNT, real* _weightsNB, real* _weightsST, real* _weightsSB);

public:
    std::vector<uint> planeNeighborNT,  planeNeighborNB, planeNeighborST, planeNeighborSB;
    std::vector<real> weightsNT, weightsNB, weightsST,  weightsSB;

protected:
    uint nPoints, nPointsRead, writingOffset;
    uint nReads=0;
};


class VTKReader : public VelocityReader
{
public:
    VTKReader(SPtr<VTKFileCollection> _fileCollection):
    fileCollection(_fileCollection)    
    {};
    void getNextVelocities(real* vx, real* vy, real* vz, real t) override;
    void getNextDistributions(real* f0, real* f1, real* f2, real* f3, real* f4, real* f5, real* f6, real* f7, real* f8, real t) override;
    void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ) override;
private:  
    uint getWriteIndex(int level, int id, int linearIdx);
    void initializeIndexVectors();

private:
    std::vector<std::vector<std::vector<uint>>> readIndices, writeIndices;
    std::vector<std::vector<size_t>> nFile;
    SPtr<VTKFileCollection> fileCollection;
};


SPtr<VelocityFileCollection> createFileCollection(std::string prefix, FileType type);
SPtr<VelocityReader> createReaderForCollection(SPtr<VelocityFileCollection> fileCollection);

#endif //VELOCITY_SETTER_H_