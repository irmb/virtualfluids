#include "BoundaryQs.h"

#include "Parameter/Parameter.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define QCOLUMNS 27

BoundaryQs::BoundaryQs(std::string path, bool isBinary)
{
	if (isBinary) 
	{
		file.open(path.c_str(), std::ios::in | std::ios::binary);
		checkFileStatus(path);
		init_Binary();
	} 
	else
	{
		file.open(path.c_str(), std::ios::in );
		checkFileStatus(path);
		init();
	}
}

BoundaryQs::BoundaryQs(std::string path, std::shared_ptr<Parameter> para, std::string str, bool isBinary)
{
	if (isBinary) file.open(path.c_str(), std::ios::in | std::ios::binary);
	else file.open(path.c_str(), std::ios::in);

	if (!file) 
		para->setObj(str, false);
	else 
	{
		para->setObj(str, true);
		if(isBinary)
			init_Binary();
		else
			init();
	}
	
}

BoundaryQs::BoundaryQs()
{

}


void BoundaryQs::checkFileStatus(std::string path)
{
	if (!file)
	{
		std::cerr << "can not open q-file: " << path << std::endl;
		exit(1);
	}
}

BoundaryQs::~BoundaryQs()
{
	file.close();
}

void BoundaryQs::init()
{
	std::vector<uint32_t> vec1D_code;
	std::string bufferString;

	file >> maxLevel;
	resizeVectors();

	for (unsigned int level = 0; level <= maxLevel; level++)
	{
		file >> levelSizes[level];
		resizeVectorsPerLevel(level, vec1D_code);
		if (levelSizes[level] == 0)
			continue;
		
		for (unsigned int elem = 0; elem < levelSizes[level]; elem++)
		{
			int columnCounter = 26;
			file >> indices[level][elem];
			file >> vec1D_code[elem];

			while (vec1D_code[elem] != 0)
			{
				if (vec1D_code[elem] % 2 == 1)
					file >> values[level][columnCounter][elem];
				vec1D_code[elem] /= 2;
				columnCounter--;
			}
			getline(file, bufferString);
		}
		vec1D_code.clear();
	}
}

void BoundaryQs::init_Binary() 
{
	std::vector<uint32_t> vec1D_code;

	std::string bufferString;
	unsigned int bufferInt;
	real bufferDouble;
	uint32_t bufferUint32_t;

	file >> maxLevel;
	resizeVectors();

	for (unsigned int level = 0; level <= maxLevel; level++)
	{
		file >> levelSizes[level];
		resizeVectorsPerLevel(level, vec1D_code);
		if (levelSizes[level] == 0)
			continue;

		for (unsigned int elem = 0; elem < levelSizes[level]; elem++)
		{
			int zaehler = 26;
			file.read((char*)&bufferInt, sizeof(int));
			indices[level][elem] = bufferInt;

			file.read((char*)&bufferUint32_t, sizeof(uint32_t));
			vec1D_code[elem] = bufferUint32_t;
			while (vec1D_code[elem] != 0)
			{
				if (vec1D_code[elem] % 2 == 1)
				{
					file.read((char*)&bufferDouble, sizeof(double));
					values[level][zaehler][elem] = bufferDouble;
				}
				vec1D_code[elem] /= 2;
				zaehler--;
			}
			getline(file, bufferString);
		}
		vec1D_code.clear();
	}
}

void BoundaryQs::resizeVectors()
{
	levelSizes.resize(maxLevel + 1);
	values.resize(maxLevel + 1);
	indices.resize(maxLevel + 1);
}

void BoundaryQs::resizeVectorsPerLevel(unsigned int level, std::vector<uint32_t> &vec1D_code)
{
	values[level].resize(QCOLUMNS);
	for (int i = 0; i < QCOLUMNS; i++)
		values[level][i].resize(levelSizes[level], -1);
	indices[level].resize(levelSizes[level]);
	vec1D_code.resize(levelSizes[level]);
}

unsigned int BoundaryQs::getSize(unsigned int level)
{
	return this->levelSizes[level];
}

unsigned int BoundaryQs::getLevel()
{
	return maxLevel;
}


void BoundaryQs::setValues(real** q27, unsigned int level) const
{
	for (int column = 0; column < values[level].size(); column++)
		for (int index = 0; index < values[level][column].size(); index++)
			q27[column][index] = values[level][column][index];
}

void BoundaryQs::setIndex(int *data, unsigned int level) const
{
	for (int index = 0; index < indices[level].size(); index++)
		data[index] = indices[level][index];
}
