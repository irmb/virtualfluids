#ifndef VELOCITY_SETTER_H_
#define VELOCITY_SETTER_H_

#include "PreCollisionInteractor.h"
#include "PrecursorWriter.h"

#include <string>
#include <vector>

class Grid;
namespace gg
{
    class BoundaryCondition;
}
class VelocityReader;

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

    void getVelocities(real* vx, real* vy, real* vz, std::vector<uint> readIndeces, std::vector<uint> writeIndices, uint offsetRead);
    bool markNANs(std::vector<uint> readIndices);
    bool inBoundingBox(real posX, real posY, real posZ){return  inXBounds(posX) && inYBounds(posY) && inZBounds(posZ); };
    bool inXBounds(real posX){ return posX<=maxX && posX>=minX; };
    bool inYBounds(real posY){ return posY<=maxY && posY>=minY; };
    bool inZBounds(real posZ){ return posZ<=maxZ && posZ>=minZ; };
    int findNeighborESB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)  , getIdxSY(posY)  , getIdxBZ(posZ)  ); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborEST(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)  , getIdxSY(posY)  , getIdxBZ(posZ)+1); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborENB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)  , getIdxSY(posY)+1, getIdxBZ(posZ)  ); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborENT(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)  , getIdxSY(posY)+1, getIdxBZ(posZ)+1); return idx>=0 && idx<nx*ny*nz ? idx : -1;};
    int findNeighborWSB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)+1, getIdxSY(posY)  , getIdxBZ(posZ)  ); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborWST(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)+1, getIdxSY(posY)  , getIdxBZ(posZ)+1); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborWNB(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)+1, getIdxSY(posY)+1, getIdxBZ(posZ)  ); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int findNeighborWNT(real posX, real posY, real posZ){ int idx = getLinearIndex(getIdxEX(posX)+1, getIdxSY(posY)+1, getIdxBZ(posZ)+1); return idx>=0 && idx<nx*ny*nz ? idx : -1; };
    int getIdxX(int linearIdx){ return linearIdx-nx*(getIdxY(linearIdx)+ny*getIdxZ(linearIdx));};
    int getIdxY(int linearIdx){ return linearIdx/nx-getIdxZ(linearIdx)*ny;};
    int getIdxZ(int linearIdx){ return linearIdx/getNumberOfPointsInXYPlane(); };
    real getX(int linearIdx){ return getIdxX(linearIdx)*deltaX+minX; };
    real getY(int linearIdx){ return getIdxY(linearIdx)*deltaY+minY; };
    real getZ(int linearIdx){ return getIdxZ(linearIdx)*deltaZ+minZ; };
    int getIdxEX(real posX){ return (posX-minX)/deltaX; };
    int getIdxSY(real posY){ return (posY-minY)/deltaY; };
    int getIdxBZ(real posZ){ return (posZ-minZ)/deltaZ; };
    int getLinearIndex(int idxX, int idxY, int idxZ){ return idxX + nx*(idxY+ny*idxZ);};
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
    std::vector<double>vxFile, vyFile, vzFile;
    bool loaded;
};


class VelocityFileCollection
{
public:
    VelocityFileCollection(std::string _prefix): 
    prefix(_prefix){};

    static VelocityFileCollection* createFileCollection(std::string prefix, std::string suffix);
    virtual VelocityReader* createReaderForCollection()=0;

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

    VelocityReader* createReaderForCollection() override;

private:
    void findFiles();
    std::string makeFileName(int level, int id, int part){ return PrecursorWriter::makeFileName(prefix, level, id, part) + ".bin." + suffix; };

public:
    static inline const std::string suffix = "vti";
    std::vector<std::vector<std::vector<VTKFile>>> files;
};


class VelocityReader
{
public:
    VelocityReader(VelocityFileCollection* _fileCollection):
    fileCollection(_fileCollection)
    { 
        this->nPoints = 0; 
        this->nPointsRead = 0;
        
    };

    virtual void getNextVelocities(real* vx, real* vy, real* vz, real t)=0;
    virtual void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ)=0;
    uint getNPoints(){return nPoints; };
    void switchLastAndNext();
    void copyIndexAndWeightToArrays();

public:
    std::vector<uint> planeNeighborNT,  planeNeighborNB, planeNeighborST, planeNeighborSB;
    std::vector<real> weightsNT, weightsNB, weightsST,  weightsSB;
    uint* planeNeighborNTH, *planeNeighborNBH, *planeNeighborSTH, *planeNeighborSBH;
    uint* planeNeighborNTD, *planeNeighborNBD, *planeNeighborSTD, *planeNeighborSBD;
    real* weightsNTH, *weightsNBH, *weightsSTH,  *weightsSBH;
    real* weightsNTD, *weightsNBD, *weightsSTD,  *weightsSBD;
    real* vxLastD, *vyLastD, *vzLastD;
    real* vxCurrentD, *vyCurrentD, *vzCurrentD;
    real* vxNextD, *vyNextD, *vzNextD;
    real* vxLastH, *vyLastH, *vzLastH;
    real* vxCurrentH, *vyCurrentH, *vzCurrentH;
    real* vxNextH, *vyNextH, *vzNextH;

protected:
    SPtr<VelocityFileCollection> fileCollection;
    uint nPoints, nPointsRead;
};


class VTKReader : public VelocityReader
{
public:
    VTKReader(VTKFileCollection* _fileCollection):
    fileCollection(_fileCollection),
    VelocityReader(_fileCollection)
    {};

    void getNextVelocities(real* vx, real* vy, real* vz, real t) override;
    void fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ) override;
private:  
    int getReadIndex(int level, int id, int linearIdx);
    int getWriteIndex(int level, int id, int linearIdx);
    void initializeIndexVectors();

private:
    std::vector<std::vector<std::vector<uint>>> readIndices, writeIndices;
    std::vector<std::vector<uint>> nFile;
    SPtr<VTKFileCollection> fileCollection;
};


class VelocitySetter : public PreCollisionInteractor
{
public:
    VelocitySetter(
        SPtr<VelocityFileCollection> _fileCollection,
        uint _NTRead):
    nTRead(_NTRead),
    fileCollection(_fileCollection)
    {    };

    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;
    void setBCArrays(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundary);

    VelocityReader* getVelocityReader(int level){ return velocityReaders[level]; };

public:
    std::vector<cudaStream_t> streams;
private:
    SPtr<VelocityFileCollection> fileCollection;
    std::vector<VelocityReader*> velocityReaders;
    std::vector<uint> nReads;
    uint nTRead;
};

#endif //VELOCITY_SETTER_H_