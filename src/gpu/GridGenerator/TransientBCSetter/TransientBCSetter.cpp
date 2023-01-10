#include "TransientBCSetter.h"
#include "GridGenerator/grid/Grid.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include <logger/Logger.h>


#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

SPtr<FileCollection> createFileCollection(std::string prefix, FileType type)
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

SPtr<TransientBCInputFileReader> createReaderForCollection(SPtr<FileCollection> fileCollection, uint readLevel)
{
    switch(fileCollection->getFileType())
    {
        case FileType::VTK:
            return std::make_shared<VTKReader>(std::static_pointer_cast<VTKFileCollection>(fileCollection), readLevel);
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

std::string readElement(std::string line)
{
    size_t elemStart = line.find("<")+1;
    // size_t elemEnd = line.find("/>", elemStart);
    size_t nameLen = line.find(" ", elemStart)-elemStart;
    return line.substr(elemStart, nameLen);
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

    getline(file, line);
    while(strcmp(readElement(line).c_str(), "DataArray")==0)
    {
        Quantity quant = Quantity();
        quant.name = readAttribute(line, "Name");
        quant.offset = std::stoi(readAttribute(line, "offset"));
        this->quantities.push_back( quant );
        getline(file, line);
    }
    getline(file, line); // </Piece
    getline(file, line); // </ImageData
    getline(file, line); // AppendedData

    int offset = int(file.tellg())+sizeof(char)+4; // skip underscore and bytesPerVal

    for(auto& quantity: this->quantities)
    {
        quantity.offset += offset;
    }

    file.close();

    this->deltaX = spacing[0];
    this->deltaY = spacing[1];
    this->deltaZ = spacing[2];

    this->nx = pieceExtent[1]-pieceExtent[0]+1;
    this->ny = pieceExtent[3]-pieceExtent[2]+1;
    this->nz = pieceExtent[5]-pieceExtent[4]+1;

    this->minX = origin[0]+this->deltaX*pieceExtent[0]; this->maxX = (this->nx-1)*this->deltaX+this->minX;
    this->minY = origin[1]+this->deltaY*pieceExtent[2]; this->maxY = (this->ny-1)*this->deltaY+this->minY;
    this->minZ = origin[2]+this->deltaZ*pieceExtent[4]; this->maxZ = (this->nz-1)*this->deltaZ+this->minZ;
    // printFileInfo();

}

bool VTKFile::markNANs(std::vector<uint> readIndices)
{
    std::ifstream buf(fileName.c_str(), std::ios::in | std::ios::binary);

    std::vector<double> tmp;
    tmp.reserve(readIndices.size());
    buf.seekg(this->quantities[0].offset);
    buf.read((char*) tmp.data(), sizeof(double)*readIndices.size());
    auto firstNAN = std::find_if(tmp.begin(), tmp.end(), [](auto it){ return isnan(it); });
    
    return firstNAN != tmp.end();
}

void VTKFile::loadFile()
{
    std::ifstream buf(this->fileName.c_str(), std::ios::in | std::ios::binary);
    for(auto& quantity: this->quantities)
    {
        quantity.values.resize(getNumberOfPoints());
        buf.seekg(quantity.offset);
        buf.read(reinterpret_cast<char*>(quantity.values.data()), this->getNumberOfPoints()*sizeof(double));
    }

    buf.close();

    this->loaded = true;
}

void VTKFile::unloadFile()
{
    for(auto& quantity : this->quantities)
    {
        std::vector<double> replacement;
        quantity.values.swap(replacement);
    }
    this->loaded = false;
}

void VTKFile::getData(real* data, uint numberOfNodes, const std::vector<uint>& readIndices, const std::vector<uint>& writeIndices, uint offsetRead, uint offsetWrite)
{
    if(!this->loaded) loadFile();

    size_t nPoints = writeIndices.size();

    for(size_t j=0; j<this->quantities.size(); j++)
    {
        real* quant = &data[j*numberOfNodes];
        for(size_t i=0; i<nPoints; i++)
        {
            quant[offsetWrite+writeIndices[i]] = this->quantities[j].values[readIndices[i]+offsetRead];
        }
    }
}

void VTKFile::printFileInfo()
{
    printf("file %s with \n nx %i ny %i nz %i \n origin %f %f %f \n spacing %f %f %f \n", 
            fileName.c_str(), nx, ny, nz, minX, minY, minZ, deltaX, deltaY, deltaZ);
    for(auto quantity: this->quantities)
    {
        printf("\t quantity %s offset %i \n", quantity.name.c_str(), quantity.offset);
    }
        
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
            {
                VF_LOG_INFO("VTKFileCollection found {} files with ID {} level {}", filesWithThisId.size(), filesOnThisLevel.size(), files.size() );
                filesOnThisLevel.push_back(filesWithThisId);
            }
            else foundLastID = true;
        }


        if(!filesOnThisLevel.empty())
            files.push_back(filesOnThisLevel);
        else 
            foundLastLevel = true;

    }

    if(files.empty())
        VF_LOG_CRITICAL("VTKFileCollection found no files!"); 
}
    
void TransientBCInputFileReader::getNeighbors(uint* neighbor0PP, uint* neighbor0PM, uint* neighbor0MP, uint* neighbor0MM)
{
    std::copy(planeNeighbor0PP.begin(), planeNeighbor0PP.end(), &neighbor0PP[writingOffset]);
    std::copy(planeNeighbor0PM.begin(), planeNeighbor0PM.end(), &neighbor0PM[writingOffset]);
    std::copy(planeNeighbor0MP.begin(), planeNeighbor0MP.end(), &neighbor0MP[writingOffset]);
    std::copy(planeNeighbor0MM.begin(), planeNeighbor0MM.end(), &neighbor0MM[writingOffset]);
}

void TransientBCInputFileReader::getWeights(real* _weights0PP, real* _weights0PM, real* _weights0MP, real* _weights0MM)
{
    std::copy(weights0PP.begin(), weights0PP.end(), &_weights0PP[writingOffset]);
    std::copy(weights0PM.begin(), weights0PM.end(), &_weights0PM[writingOffset]);
    std::copy(weights0MP.begin(), weights0MP.end(), &_weights0MP[writingOffset]);
    std::copy(weights0MM.begin(), weights0MM.end(), &_weights0MM[writingOffset]);
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

    this->weights0PP.reserve(this->nPoints);
    this->weights0PM.reserve(this->nPoints);
    this->weights0MP.reserve(this->nPoints);
    this->weights0MM.reserve(this->nPoints);

    this->planeNeighbor0PP.reserve(this->nPoints);
    this->planeNeighbor0PM.reserve(this->nPoints);
    this->planeNeighbor0MP.reserve(this->nPoints);
    this->planeNeighbor0MM.reserve(this->nPoints);

    for(uint i=0; i<nPoints; i++)
    {

        real posY = coordsY[i];
        real posZ = coordsZ[i];
        bool found0PP = false, found0PM = false, found0MP = false, found0MM = false, foundAll = false;

        uint level = this->readLevel;

        for(int fileId=0; fileId<(int)this->fileCollection->files[level].size(); fileId++)
        {
            VTKFile &file = this->fileCollection->files[level][fileId][0];
            if(!file.inBoundingBox(posY, posZ, 0.0f)) continue;

            // y in simulation is x in precursor/file, z in simulation is y in precursor/file 
            // simulation -> file: N -> E, S -> W, T -> N, B -> S
            int idx = file.findNeighborMMM(posY, posZ, 0.f);                            //!> index of nearest WSB neighbor on precursor file
            
            if(idx!=-1)
            {
                // Filter for exact matches
                if(abs(posY-file.getX(idx)) < max_diff && abs(posZ-file.getY(idx)) < max_diff) 
                {
                    this->weights0PP.emplace_back(1e6f);
                    this->weights0PM.emplace_back(0.f);
                    this->weights0MP.emplace_back(0.f);
                    this->weights0MM.emplace_back(0.f);
                    uint writeIdx = this->getWriteIndex(level, fileId, idx);            //!> writeIdx: index on host/device array where precursor value will be written to after loading from file
                    this->planeNeighbor0PP.push_back(writeIdx);                          //!> neighbor lists mapping where BC kernel should read from on host/device array
                    this->planeNeighbor0PM.push_back(writeIdx);
                    this->planeNeighbor0MP.push_back(writeIdx);
                    this->planeNeighbor0MM.push_back(writeIdx);
                    found0PP = true; found0PM = true; found0MM = true; found0MP = true;
                } 
                else
                {
                    perfect_match = false;
                }

                if(!found0MM)
                {
                    found0MM = true;
                    real dy = file.getX(idx)-posY;
                    real dz = file.getY(idx)-posZ;
                    this->weights0MM.emplace_back(1.f/(dy*dy+dz*dz+eps));
                    this->planeNeighbor0MM.emplace_back(getWriteIndex(level, fileId, idx));
                }
                
            } 
            
            if(!found0PP) //NT in simulation is EN in precursor
            {
                int idx = file.findNeighborPPM(posY, posZ, 0.f);
                if(idx!=-1)
                {
                    found0PP = true;
                    real dy = file.getX(idx)-posY;
                    real dz = file.getY(idx)-posZ;
                    this->weights0PP.emplace_back(1.f/(dy*dy+dz*dz+eps));
                    this->planeNeighbor0PP.emplace_back(getWriteIndex(level, fileId, idx));
                }
            }

            if(!found0PM) //NB in simulation is ES in precursor
            {
                int idx = file.findNeighborPMM(posY, posZ, 0.f);
                if(idx!=-1)
                {
                    found0PM = true;
                    real dy = file.getX(idx)-posY;
                    real dz = file.getY(idx)-posZ;
                    this->weights0PM.emplace_back(1.f/(dy*dy+dz*dz+eps));
                    this->planeNeighbor0PP.emplace_back(getWriteIndex(level, fileId, idx));
                }
            }

            if(!found0MP) //ST in simulation is WN in precursor
            {
                int idx = file.findNeighborMPM(posY, posZ, 0.f);
                if(idx!=-1)
                {
                    found0MP = true;
                    real dy = file.getX(idx)-posY;
                    real dz = file.getY(idx)-posZ;
                    this->weights0MP.emplace_back(1.f/(dy*dy+dz*dz+eps));
                    this->planeNeighbor0MP.emplace_back(getWriteIndex(level, fileId, idx));
                }
            }

            foundAll = found0PP && found0PM && found0MP && found0MM;

            if(foundAll) break;
        }

        if(!foundAll)
        {
            VF_LOG_CRITICAL("Found no matching precursor neighbors for grid point at y={}, z={} \n", posY, posZ);
            throw std::runtime_error("VTKReader::fillArrays(): Did not find neighbors in the FileCollection for all points");
        }
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
        this->writeIndices[level][id].push_back(this->nPointsRead);     //!> index on host/device array where value from file will be written to
        this->readIndices[level][id].push_back(linearIndex);            //!> index in file that will be read from 
        this->nPointsRead++;
    }
    return idx;
}


void VTKReader::getNextData(real* data, uint numberOfNodes, real time)
{
    // for(size_t level=0; level<this->fileCollection->files.size(); level++)
    // {
        uint level = this->readLevel;
        for(size_t id=0; id<this->fileCollection->files[level].size(); id++)
        {
            size_t numberOfFiles = this->nFile[level][id];


            if(!this->fileCollection->files[level][id][numberOfFiles].inZBounds(time))
            {
                numberOfFiles++;

                VF_LOG_INFO("PrecursorBC on level {}: switching to file no. {}\n", level, numberOfFiles);
                if(numberOfFiles == this->fileCollection->files[level][id].size())
                    throw std::runtime_error("Not enough Precursor Files to read");

                this->fileCollection->files[level][id][numberOfFiles-1].unloadFile();
                if(numberOfFiles+1<this->fileCollection->files[level][id].size())
                {
                    VTKFile* nextFile = &this->fileCollection->files[level][id][numberOfFiles+1];
                    if(! nextFile->isLoaded())
                    {
                        read.wait();
                        read = std::async(std::launch::async, [](VTKFile* file){ file->loadFile(); }, &this->fileCollection->files[level][id][numberOfFiles+1]);
                    }
                }
            }
        

            VTKFile* file = &this->fileCollection->files[level][id][numberOfFiles];

            int off = file->getClosestIdxZ(time)*file->getNumberOfPointsInXYPlane();
            file->getData(data, numberOfNodes, this->readIndices[level][id], this->writeIndices[level][id], off, this->writingOffset);
            this->nFile[level][id] = numberOfFiles;
        }
    // }
}