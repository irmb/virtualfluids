#include "CoordNeighborGeoV.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

CoordNeighborGeoV::CoordNeighborGeoV()
{

}

CoordNeighborGeoV::CoordNeighborGeoV(std::string path, bool binaer, bool coord)
{
	file.open(path.c_str(), std::ios::in | std::ios::binary);
	if (!file) 
	{
		std::cerr << "error: can not open file CoordNeighborGeo: " << path << std::endl;
		exit(1);
	}
	if(binaer)
		init_Binary(coord);
	else
		init(coord);	
}

CoordNeighborGeoV::~CoordNeighborGeoV()
{
	file.close();
}

void CoordNeighborGeoV::init(bool isCoord)
{
	this->readLevel();
	resizeVectors();

	for (unsigned int i = 0; i <= maxLevel; i++)
	{
		readLevelSize(i);
		if (isCoord)
		{
			coordinates[i].resize(levelSizes[i] + 1);
			for (unsigned int j = 0; j <= levelSizes[i]; j++)
				file >> coordinates[i][j];
		}
		else
		{
			neighbors[i].resize(levelSizes[i] + 1);
			for (unsigned int j = 0; j <= levelSizes[i]; j++)
				file >> neighbors[i][j];
		}
	}
}


void CoordNeighborGeoV::resizeVectors()
{
	levelSizes.resize(maxLevel + 1);
	coordinates.resize(maxLevel + 1);
	neighbors.resize(maxLevel + 1);
}

void CoordNeighborGeoV::readLevel()
{
	file >> this->maxLevel;
}

void CoordNeighborGeoV::init_Binary(bool isCoord) 
{
	this->readLevel();
	this->resizeVectors();

	readLevelSize(0);

	for(unsigned int level = 0; level <= maxLevel; level++) 
	{
        //readLevelSize(level);
		if(isCoord) 
			readCoordinates(level);
		else 
			readNeighbors(level);

		if(level == maxLevel) break;

		skipSpace();
		readLevelSize(level+1);
	}
}

void CoordNeighborGeoV::readLevelSize(unsigned int level)
{
	file >> levelSizes[level];
}

void CoordNeighborGeoV::readNeighbors(unsigned int level)
{
	unsigned int bufferInt;
	neighbors[level].resize(levelSizes[level] + 1);
	file >> neighbors[level][0];
	skipSpace();

	for (unsigned int j = 0; j < levelSizes[level]; j++)
	{
		file.read((char*)&bufferInt, sizeof(unsigned int));
		neighbors[level][j + 1] = bufferInt;
	}
}

void CoordNeighborGeoV::readCoordinates(unsigned int level)
{
	double bufferDouble;
	coordinates[level].resize(levelSizes[level] + 1);
	file >> coordinates[level][0];
	skipSpace();

	for (unsigned int j = 0; j < levelSizes[level]; j++)
	{
		file.read((char*)&bufferDouble, sizeof(double));
		coordinates[level][j + 1] = (real)bufferDouble;
	}
}

unsigned int CoordNeighborGeoV::getLevel() 
{
	return maxLevel;
}

unsigned int CoordNeighborGeoV::getSize(unsigned int level) 
{
	return this->levelSizes[level];
}


std::vector<unsigned int> CoordNeighborGeoV::getVec(unsigned int level) 
{
	return this->neighbors[level];
}

void CoordNeighborGeoV::setVec(unsigned int level, std::vector<unsigned int> vec) {
	this->neighbors[level]=vec;
	//for (int i=0; i<=2200; i++) {
		//std::cout <<"Test im Setter: "<< i <<": " << vec[i] << std::endl;
	//}
}


void CoordNeighborGeoV::initalCoords(real *data, unsigned int level) const
{
	for (std::size_t index = 0; index < coordinates[level].size(); index++)
		data[index] = coordinates[level][index];
}

void CoordNeighborGeoV::initalNeighbors(unsigned int *data, unsigned int level) const
{
	for (std::size_t index = 0; index < neighbors[level].size(); index++)
		data[index] = neighbors[level][index];
}

void CoordNeighborGeoV::skipSpace()
{
	char c;
	file.get(c);
}
