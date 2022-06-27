#include "VelocitySetter.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GeometryUtils.h"
#include "Parameter/Parameter.h"
#include "GridGenerator/grid/Grid.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include <cuda/CudaGrid.h>
#include <logger/Logger.h>


#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
//TODO test if larger memcopy is faster
__global__ void interpolateVelocities(  uint nVelPoints,
                                        real* vxLast, real* vyLast, real* vzLast,
                                        real* vxNext, real* vyNext, real* vzNext,
                                        real* vx, real* vy, real* vz,
                                        uint* planeNeighborNT, uint* planeNeighborNB,
                                        uint* planeNeighborST, uint* planeNeighborSB,
                                        real* weightsNT, real* weightsNB,
                                        real* weightsST, real* weightsSB,
                                        real tRatio, real velocityRatio)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;
    if (node>nVelPoints) return;

    real vxLastInterpd, vyLastInterpd, vzLastInterpd; 
    real vxNextInterpd, vyNextInterpd, vzNextInterpd; 

    uint kNT = planeNeighborNT[node];
    real dNT = weightsNT[node];
    if(dNT < 1e6)
    {
        uint kNB = planeNeighborNB[node];
        uint kST = planeNeighborST[node];
        uint kSB = planeNeighborSB[node];

        real dNB = weightsNB[node];
        real dST = weightsST[node];
        real dSB = weightsSB[node];

        real invWeightSum = 1/(dNT+dNB+dST+dSB);

        vxLastInterpd = (vxLast[kNT]*dNT + vxLast[kNB]*dNB + vxLast[kST]*dST + vxLast[kSB]*dSB)*invWeightSum;
        vyLastInterpd = (vyLast[kNT]*dNT + vyLast[kNB]*dNB + vyLast[kST]*dST + vyLast[kSB]*dSB)*invWeightSum;
        vzLastInterpd = (vzLast[kNT]*dNT + vzLast[kNB]*dNB + vzLast[kST]*dST + vzLast[kSB]*dSB)*invWeightSum;

        vxNextInterpd = (vxNext[kNT]*dNT + vxNext[kNB]*dNB + vxNext[kST]*dST + vxNext[kSB]*dSB)*invWeightSum;
        vyNextInterpd = (vyNext[kNT]*dNT + vyNext[kNB]*dNB + vyNext[kST]*dST + vyNext[kSB]*dSB)*invWeightSum;
        vzNextInterpd = (vzNext[kNT]*dNT + vzNext[kNB]*dNB + vzNext[kST]*dST + vzNext[kSB]*dSB)*invWeightSum;
    }
    else
    {
        vxLastInterpd = vxLast[kNT];
        vyLastInterpd = vyLast[kNT];
        vzLastInterpd = vzLast[kNT];

        vxNextInterpd = vxNext[kNT];
        vyNextInterpd = vyNext[kNT];
        vzNextInterpd = vzNext[kNT];
    }
    // if(node==100)
        // printf("last u %f v %f next u %f v %f\n", vxLastInterpd, vyLastInterpd, vxNextInterpd, vyNextInterpd);
    vx[node] = ((1.f-tRatio)*vxLastInterpd + tRatio*vxNextInterpd)/velocityRatio;
    vy[node] = ((1.f-tRatio)*vyLastInterpd + tRatio*vyNextInterpd)/velocityRatio; 
    vz[node] = ((1.f-tRatio)*vzLastInterpd + tRatio*vzNextInterpd)/velocityRatio;
}

template<typename T>
std::vector<T> readStringToVector(std::string s)
{
    std::vector<T> out;
    std::stringstream input(s);
    float num;
    while(input >> num)
    {
        out.push_back(num);
    }
    return out;
}

std::string readAttribute(std::string line, std::string attributeName)
{
    int attributeStart = line.find(attributeName)+attributeName.size()+2; // add to for '="'
    int attributeLen = line.find("\"", attributeStart)-attributeStart;
    return line.substr(attributeStart, attributeLen);
}

SPtr<VelocityFileCollection> VelocityFileCollection::createFileCollection(std::string prefix, std::string suffix)
{
    if(strcmp(suffix.c_str(), VTKFileCollection::suffix.c_str())==0)
        return std::make_shared<VTKFileCollection>(prefix);
    else return nullptr;
}



void VTKFile::readHeader()
{
    //TODO make this more flexible
    std::ifstream file(this->fileName);

    std::string line;

    getline(file, line); // VTKFile
    if(line[1]=='?') getline(file, line); // ignore first line if xml version

    getline(file, line); // ImageData
    std::vector<int> wholeExtent = readStringToVector<int>(readAttribute(line, "WholeExtent"));
    std::vector<float> origin = readStringToVector<float>(readAttribute(line, "Origin"));
    std::vector<float> spacing = readStringToVector<float>(readAttribute(line, "Spacing"));

    getline(file, line); // Piece 
    std::vector<int> pieceExtent = readStringToVector<int>(readAttribute(line, "Extent"));
    getline(file, line); // PointData

    getline(file, line); // DataArray vx
    if(strcmp(readAttribute(line, "Name").c_str(), "vx")==0) // just to make sure
        this->offsetVx = std::stoi(readAttribute(line, "offset"));

    getline(file, line); // DataArray vy
    if(strcmp(readAttribute(line, "Name").c_str(), "vy")==0) // just to make sure
        this->offsetVy = std::stoi(readAttribute(line, "offset"));

    getline(file, line); // DataArray vz
    if(strcmp(readAttribute(line, "Name").c_str(), "vz")==0) // just to make sure
        this->offsetVz = std::stoi(readAttribute(line, "offset"));

    this->deltaX = spacing[0];
    this->deltaY = spacing[1];
    this->deltaZ = spacing[2];

    this->nx = pieceExtent[1]-pieceExtent[0]+1;
    this->ny = pieceExtent[3]-pieceExtent[2]+1;
    this->nz = pieceExtent[5]-pieceExtent[4]+1;

    this->minX = origin[0]+this->deltaX*pieceExtent[0]; this->maxX = this->nx*this->deltaX+this->minX;
    this->minY = origin[1]+this->deltaY*pieceExtent[2]; this->maxY = this->ny*this->deltaY+this->minY;
    this->minZ = origin[2]+this->deltaZ*pieceExtent[4]; this->maxZ = this->nz*this->deltaZ+this->minZ;

    getline(file, line); // </ PointData
    getline(file, line); // </Piece
    getline(file, line); // </ImageData
    getline(file, line); // AppendedData

    int offset = int(file.tellg())+sizeof(char)+4; // skip underscore and bytesPerVal

    this->offsetVx += offset;
    this->offsetVy += offset;
    this->offsetVz += offset;

    file.close();

}

bool VTKFile::markNANs(std::vector<uint> readIndices)
{
    std::ifstream buf(fileName.c_str(), std::ios::in | std::ios::binary);

    int arraySize = sizeof(double)*readIndices.size();
    std::vector<double> tmp;
    tmp.reserve(readIndices.size());
    buf.seekg(this->offsetVx);
    buf.read((char*) tmp.data(), sizeof(double)*readIndices.size());
    auto firstNAN = std::find_if(tmp.begin(), tmp.end(), [](auto it){ return isnan(it); });
    
    return firstNAN != tmp.end();
}

void VTKFile::loadFile()
{
    std::ifstream buf(this->fileName.c_str(), std::ios::in | std::ios::binary);
    this->vxFile.reserve(getNumberOfPoints()); 
    buf.seekg(this->offsetVx);
    buf.read(reinterpret_cast<char*>(this->vxFile.data()), this->getNumberOfPoints()*sizeof(double));

    this->vyFile.reserve(getNumberOfPoints()); 
    buf.seekg(this->offsetVy);
    buf.read((char*) this->vyFile.data(), this->getNumberOfPoints()*sizeof(double));

    this->vzFile.reserve(this->getNumberOfPoints());
    buf.seekg(this->offsetVz);
    buf.read((char*) this->vzFile.data(), this->getNumberOfPoints()*sizeof(double));

    buf.close();
    this->loaded = true;
}

void VTKFile::unloadFile()
{
    std::vector<double> replaceX, replaceY, replaceZ;
    this->vxFile.swap(replaceX);
    this->vyFile.swap(replaceY);
    this->vzFile.swap(replaceZ);
    this->loaded = false;
}

void VTKFile::getVelocities(real* vx, real* vy, real* vz, std::vector<uint> readIndeces, std::vector<uint> writeIndices, uint offsetRead, uint offsetWrite)
{
    if(!this->loaded) loadFile();

    uint nPoints = writeIndices.size();

    for(int i=0; i<nPoints; i++)
    {   
        vx[offsetWrite+writeIndices[i]] = this->vxFile[readIndeces[i]+offsetRead];
        vy[offsetWrite+writeIndices[i]] = this->vyFile[readIndeces[i]+offsetRead];
        vz[offsetWrite+writeIndices[i]] = this->vzFile[readIndeces[i]+offsetRead];
    }
}

void VTKFile::printFileInfo()
{
    printf("file %s with \n nx %i ny %i nz %i \n origin %f %f %f \n spacing %f %f %f \n offset %i %i %i \n", 
            fileName.c_str(), nx, ny, nz, minX, minY, minZ, deltaX, deltaY, deltaZ, offsetVx, offsetVy, offsetVz);
        
}


void VTKFileCollection::findFiles()
{
    bool foundLastLevel = false;

    while(!foundLastLevel)
    {
        bool foundLastID = false;
        std::vector<std::vector<VTKFile>> filesOnThisLevel;
        while(!foundLastID)
        {
            bool foundLastPart = false;
            std::vector<VTKFile> filesWithThisId;
            while (!foundLastPart)
            {
                std::string fname = makeFileName(files.size(), filesOnThisLevel.size(), filesWithThisId.size());
                std::ifstream f(fname);
                if(f.good())
                    filesWithThisId.emplace_back(fname);
                else
                    foundLastPart = true;
                
            }
            if(!filesWithThisId.empty())
                filesOnThisLevel.push_back(filesWithThisId);
            else foundLastID = true;
        }
        if(!filesOnThisLevel.empty())
            files.push_back(filesOnThisLevel);
        else foundLastLevel = true;
    }
}


SPtr<VelocityReader> VTKFileCollection::createReaderForCollection()
{
    return std::make_shared<VTKReader>(getSelf()); 
}


void VelocityReader::switchLastAndNext()
{
    real *tmp = vxLastD;
    vxLastD = vxCurrentD;
    vxCurrentD = vxNextD;
    vxNextD = tmp;

    tmp = vyLastD;
    vyLastD = vyCurrentD;
    vyCurrentD = vyNextD;
    vyNextD = tmp;

    tmp = vzLastD;
    vzLastD = vzCurrentD;
    vzCurrentD = vzNextD;
    vzNextD = tmp;
}

void VelocityReader::getNeighbors(uint* neighborNT, uint* neighborNB, uint* neighborST, uint* neighborSB)
{
    std::copy(planeNeighborNT.begin(), planeNeighborNT.end(), &neighborNT[writingOffset]);
    std::copy(planeNeighborNB.begin(), planeNeighborNB.end(), &neighborNB[writingOffset]);
    std::copy(planeNeighborST.begin(), planeNeighborST.end(), &neighborST[writingOffset]);
    std::copy(planeNeighborSB.begin(), planeNeighborSB.end(), &neighborSB[writingOffset]);
}

void VelocityReader::getWeights(real* _weightsNT, real* _weightsNB, real* _weightsST, real* _weightsSB)
{
    std::copy(weightsNT.begin(), weightsNT.end(), &_weightsNT[writingOffset]);
    std::copy(weightsNB.begin(), weightsNB.end(), &_weightsNB[writingOffset]);
    std::copy(weightsST.begin(), weightsST.end(), &_weightsST[writingOffset]);
    std::copy(weightsSB.begin(), weightsSB.end(), &_weightsSB[writingOffset]);
}

void VTKReader::initializeIndexVectors()
{
    this->readIndices.resize(this->fileCollection->files.size());
    this->writeIndices.resize(this->fileCollection->files.size());
    this->nFile.resize(this->fileCollection->files.size());
    for(int lev=0; lev<this->fileCollection->files.size(); lev++)
    {
        this->readIndices[lev].resize(this->fileCollection->files[lev].size());
        this->writeIndices[lev].resize(this->fileCollection->files[lev].size());
        this->nFile[lev].resize(this->fileCollection->files[lev].size());
    }
}

void VTKReader::fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ)
{
    this->nPoints = coordsY.size();
    this->initializeIndexVectors();
    real max_diff = 1e-4; // maximum distance between point on grid and precursor plane to count as exact match
    bool perfect_match = true;

    //TODO reserve vectors;

    for(uint i=0; i<nPoints; i++)
    {

        real posY = coordsY[i];
        real posZ = coordsZ[i];

        // printf("y %f z %f i %d npoints %u\n", posY, posZ, i, nPoints);
        bool foundNT = false, foundNB = false, foundST = false, foundSB = false, foundAll = false;

        for(int level=this->fileCollection->files.size()-1; level>=0; level--) // go backwards to find finest nodes first
        {
            for(int id=0; id<this->fileCollection->files[level].size(); id++)
            {
                VTKFile file = this->fileCollection->files[level][id][0];

                // Filter for exact matches
                int idx = file.findNeighborESB(posY, posZ, 0.f);
                if(idx!=-1)
                {
                    if(abs(posY<file.getX(idx))<max_diff && abs(posZ<file.getY(idx))< max_diff)
                    {
                        this->weightsNT.emplace_back(1e6f);
                        this->weightsNB.emplace_back(0.f);
                        this->weightsST.emplace_back(0.f);
                        this->weightsSB.emplace_back(0.f);
                        uint writeIdx = this->getWriteIndex(level, id, idx);
                        this->planeNeighborNT.push_back(writeIdx);
                        this->planeNeighborNB.push_back(writeIdx);
                        this->planeNeighborST.push_back(writeIdx);
                        this->planeNeighborSB.push_back(writeIdx);
                        foundNT = true; foundNB = true; foundSB = true; foundST = true;
                    }
                } else perfect_match = false;



                if(!foundNT)
                {
                    int idx = file.findNeighborENT(posY, posZ, 0.f); // VTKFile is x,y,z but VTKCollection is y,z,t
                    if(idx!=-1)
                    {
                        foundNT = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsNT.emplace_back(1.f/(dy*dy+dz*dz+1e-7));
                        this->planeNeighborNT.emplace_back(getWriteIndex(level, id, idx));
                    }

                }
                if(!foundNB)
                {
                    int idx = file.findNeighborEST(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundNB = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsNB.emplace_back(1.f/(dy*dy+dz*dz+1e-7));
                        this->planeNeighborNT.emplace_back(getWriteIndex(level, id, idx));
                    }
                }
                if(!foundST)
                {
                    int idx = file.findNeighborWNT(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundST = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsST.emplace_back(1.f/(dy*dy+dz*dz+1e-7));
                        this->planeNeighborST.emplace_back(getWriteIndex(level, id, idx));
                    }
                }
                if(!foundSB)
                {
                    int idx = file.findNeighborWST(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundSB = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsSB.emplace_back(1.f/(dy*dy+dz*dz+1e-7));
                        this->planeNeighborSB.emplace_back(getWriteIndex(level, id, idx));
                    }
                }

                foundAll = foundNT && foundNB && foundST && foundSB;

                if(foundAll) break;
            }
            if(foundAll) break;
        }
        if(!foundAll)
            throw std::runtime_error("Did not find neighbors in the VelocityFileCollection for all points");
    }

    if(perfect_match)
        printf("was a perfect match \n");

    for(int level=0; level<this->fileCollection->files.size(); level++){
        for(int id=0; id<this->fileCollection->files[level].size(); id++){
            if(this->fileCollection->files[level][id][0].markNANs(this->readIndices[level][id]))
                throw std::runtime_error("Found a NAN in the precursor where a velocity is needed");
    }}
}

int VTKReader::getWriteIndex(int level, int id, int linearIndex)
{
    auto it = std::find(this->writeIndices[level][id].begin(), this->writeIndices[level][id].end(), linearIndex);
    int idx = it-this->writeIndices[level][id].begin();
    if(it==this->writeIndices[level][id].end())
    {
        this->writeIndices[level][id].push_back(this->nPointsRead);
        this->readIndices[level][id].push_back(linearIndex);
        this->nPointsRead++;
    }
    return idx;
}

void VTKReader::getNextVelocities(real* vx, real* vy, real* vz, real t)
{
    for(int level=0; level<this->fileCollection->files.size(); level++)
    {
        for(int id=0; id<this->fileCollection->files[level].size(); id++)
        {
            int nF = this->nFile[level][id];
            if(!this->fileCollection->files[level][id][nF].inZBounds(t))
            {
                this->fileCollection->files[level][id][nF].unloadFile();
                nF++;
            }
        
            if(nF == this->fileCollection->files[level][id].size())
                throw std::runtime_error("Not enough Precursor Files to read");

            auto file = &this->fileCollection->files[level][id][nF];

            int off = file->getIdxBZ(t)*file->getNumberOfPointsInXYPlane();
            file->getVelocities(vx, vy, vz, this->readIndices[level][id], this->writeIndices[level][id], off, this->writingOffset);
            this->nFile[level][id] = nF;
        }
    }
}