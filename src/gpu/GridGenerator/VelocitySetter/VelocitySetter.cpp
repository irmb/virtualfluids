#include "VelocitySetter.h"
#include "GridGenerator/grid/Grid.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include <logger/Logger.h>


#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

SPtr<VelocityFileCollection> createFileCollection(std::string prefix, FileType type)
{
    switch(type)
    {
        case FileType::VTK:
            return std::make_shared<VTKFileCollection>(prefix);
            break;
        default:
            return nullptr;
    }
}

SPtr<VelocityReader> createReaderForCollection(SPtr<VelocityFileCollection> fileCollection)
{
    switch(fileCollection->getFileType())
    {
        case FileType::VTK:
            return std::make_shared<VTKReader>(std::static_pointer_cast<VTKFileCollection>(fileCollection));
            break;
        default:
            return nullptr;
    }
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
    size_t attributeStart = line.find(attributeName)+attributeName.size() + 2; // add 2 for '="'
    size_t attributeLen = line.find("\"", attributeStart)-attributeStart;
    return line.substr(attributeStart, attributeLen);
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

    printFileInfo();

}

bool VTKFile::markNANs(std::vector<uint> readIndices)
{
    std::ifstream buf(fileName.c_str(), std::ios::in | std::ios::binary);

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
    this->vxFile.resize(getNumberOfPoints()); 
    buf.seekg(this->offsetVx);
    buf.read(reinterpret_cast<char*>(this->vxFile.data()), this->getNumberOfPoints()*sizeof(double));

    this->vyFile.resize(getNumberOfPoints()); 
    buf.seekg(this->offsetVy);
    buf.read((char*) this->vyFile.data(), this->getNumberOfPoints()*sizeof(double));

    this->vzFile.resize(this->getNumberOfPoints());
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

    size_t nPoints = writeIndices.size();

    for(size_t i=0; i<nPoints; i++)
    {   
        // if(i<10)printf("point read %d point write %d idx read %d idx write %d\n", readIndeces[i], writeIndices[i], readIndeces[i]+offsetRead, offsetWrite+writeIndices[i]);
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
                std::string fname = makeFileName((int)files.size(), (int)filesOnThisLevel.size(), (int)filesWithThisId.size());
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
    for(size_t lev=0; lev<this->fileCollection->files.size(); lev++)
    {
        this->readIndices[lev].resize(this->fileCollection->files[lev].size());
        this->writeIndices[lev].resize(this->fileCollection->files[lev].size());
        this->nFile[lev].resize(this->fileCollection->files[lev].size());
    }
}

void VTKReader::fillArrays(std::vector<real>& coordsY, std::vector<real>& coordsZ)
{
    this->nPoints = (uint)coordsY.size();
    this->initializeIndexVectors();
    real max_diff = 1e-4; // maximum distance between point on grid and precursor plane to count as exact match
    real eps = 1e-7; // small number to avoid division by zero
    bool perfect_match = true;

    this->weightsNT.reserve(this->nPoints);
    this->weightsNB.reserve(this->nPoints);
    this->weightsST.reserve(this->nPoints);
    this->weightsSB.reserve(this->nPoints);

    this->planeNeighborNT.reserve(this->nPoints);
    this->planeNeighborNB.reserve(this->nPoints);
    this->planeNeighborST.reserve(this->nPoints);
    this->planeNeighborSB.reserve(this->nPoints);

    for(uint i=0; i<nPoints; i++)
    {

        real posY = coordsY[i];
        real posZ = coordsZ[i];
        bool foundNT = false, foundNB = false, foundST = false, foundSB = false, foundAll = false;

        for(int level= static_cast<int>(this->fileCollection->files.size())-1; level>=0; level--) // go backwards to find finest nodes first
        {
            for(size_t fileId=0; fileId<this->fileCollection->files[level].size(); fileId++)
            {
                VTKFile file = this->fileCollection->files[level][fileId][0];
                if(!file.inBoundingBox(posY, posZ, 0.0f)) continue;
                // y in simulation is x in precursor/file, z in simulation is y in precursor/file 
                // simulation -> file: N -> E, S -> W, T -> N, B -> S
                int idx = file.findNeighborWSB(posY, posZ, 0.f);
                if(idx!=-1)
                {
                    // Filter for exact matches
                    if(abs(posY-file.getX(idx)) < max_diff && abs(posZ-file.getY(idx)) < max_diff) 
                    {
                        this->weightsNT.emplace_back(1e6f);
                        this->weightsNB.emplace_back(0.f);
                        this->weightsST.emplace_back(0.f);
                        this->weightsSB.emplace_back(0.f);
                        uint writeIdx = this->getWriteIndex(level, fileId, idx);
                        this->planeNeighborNT.push_back(writeIdx);
                        this->planeNeighborNB.push_back(writeIdx);
                        this->planeNeighborST.push_back(writeIdx);
                        this->planeNeighborSB.push_back(writeIdx);
                        foundNT = true; foundNB = true; foundSB = true; foundST = true;
                    } 
                    else
                    {
                        perfect_match = false;
                    }

                    if(!foundSB)
                    {
                        foundSB = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsSB.emplace_back(1.f/(dy*dy+dz*dz+eps));
                        this->planeNeighborSB.emplace_back(getWriteIndex(level, fileId, idx));
                    }
                    
                } 

                if(!foundNT) //NT in simulation is EN in precursor
                {
                    int idx = file.findNeighborENB(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundNT = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsNT.emplace_back(1.f/(dy*dy+dz*dz+eps));
                        this->planeNeighborNT.emplace_back(getWriteIndex(level, fileId, idx));
                    }
                }

                if(!foundNB) //NB in simulation is ES in precursor
                {
                    int idx = file.findNeighborESB(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundNB = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsNB.emplace_back(1.f/(dy*dy+dz*dz+eps));
                        this->planeNeighborNT.emplace_back(getWriteIndex(level, fileId, idx));
                    }
                }

                if(!foundST) //ST in simulation is WN in precursor
                {
                    int idx = file.findNeighborWNB(posY, posZ, 0.f);
                    if(idx!=-1)
                    {
                        foundST = true;
                        real dy = file.getX(idx)-posY;
                        real dz = file.getY(idx)-posZ;
                        this->weightsST.emplace_back(1.f/(dy*dy+dz*dz+eps));
                        this->planeNeighborST.emplace_back(getWriteIndex(level, fileId, idx));
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
        printf("Precursor was a perfect match \n");


    for(size_t level=0; level<this->fileCollection->files.size(); level++){
        for(size_t id=0; id<this->fileCollection->files[level].size(); id++){
            if(this->fileCollection->files[level][id][0].markNANs(this->readIndices[level][id]))
                throw std::runtime_error("Found a NAN in the precursor where a velocity is needed");
    }}
}

uint VTKReader::getWriteIndex(int level, int id, int linearIndex)
{
    auto it = std::find(this->writeIndices[level][id].begin(), this->writeIndices[level][id].end(), linearIndex);
    uint idx = it-this->writeIndices[level][id].begin();
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
    for(size_t level=0; level<this->fileCollection->files.size(); level++)
    {
        for(size_t id=0; id<this->fileCollection->files[level].size(); id++)
        {
            size_t nF = this->nFile[level][id];
            if(!this->fileCollection->files[level][id][nF].inZBounds(t))
            {
                this->fileCollection->files[level][id][nF].unloadFile();
                nF++;
                printf("switching to file %ld\n", nF);
            }
        
            if(nF == this->fileCollection->files[level][id].size())
                throw std::runtime_error("Not enough Precursor Files to read");

            VTKFile* file = &this->fileCollection->files[level][id][nF];

            int off = (file->getClosestIdxZ(t))*file->getNumberOfPointsInXYPlane();
            file->getVelocities(vx, vy, vz, this->readIndices[level][id], this->writeIndices[level][id], off, this->writingOffset);
            this->nFile[level][id] = nF;
        }
    }
}