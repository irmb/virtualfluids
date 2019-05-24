#include "StreetPointFinder.h"

#include "Core/Logger/Logger.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "grid/Grid.h"
#include "grid/NodeValues.h"

Street::Street(real xStartCell, real yStartCell, real xEndCell, real yEndCell, real dx)
{
	real length = std::sqrt((xEndCell - xStartCell)*(xEndCell - xStartCell)
		+ (yEndCell - yStartCell)*(yEndCell - yStartCell));

	this->numberOfCells = std::floor(length / dx);

	real realLength = dx * (this->numberOfCells);

	real vectorX = (xEndCell - xStartCell) / length;
	real vectorY = (yEndCell - yStartCell) / length;

	this->xStart = xStartCell - 0.5 * (realLength - length) * vectorX + 0.5 * dx  * vectorX;
	this->yStart = yStartCell - 0.5 * (realLength - length) * vectorY + 0.5 * dx  * vectorY;

	this->xEnd = xEndCell + 0.5 * (realLength - length) * vectorX - 0.5 * dx  * vectorX;
	this->yEnd = yEndCell + 0.5 * (realLength - length) * vectorY - 0.5 * dx  * vectorY;

	//this->xStart = xStart + dx * (xEnd - xStart) / length;
	//this->yStart = yStart + dx * (yEnd - yStart) / length;
	//
	//this->xEnd   = xEnd   - dx * (xEnd - xStart) / length;
	//this->yEnd   = yEnd   - dx * (yEnd - yStart) / length;

	//this->numberOfCells = std::lround( length / dx ) + 1;
}

real Street::getCoordinateX(int cellIndex)
{
	return xStart + real(cellIndex) / real(numberOfCells - 1) * (xEnd - xStart);
}

real Street::getCoordinateY(int cellIndex)
{
	return yStart + real(cellIndex) / real(numberOfCells - 1) * (yEnd - yStart);
}

real Street::getVectorX()
{
	real vecX = this->xEnd - this->xStart;
	real vecY = this->yEnd - this->yStart;

	real length = sqrt(vecX*vecX + vecY*vecY);

	return vecX / length;
}

real Street::getVectorY()
{
	real vecX = this->xEnd - this->xStart;
	real vecY = this->yEnd - this->yStart;

	real length = sqrt(vecX*vecX + vecY*vecY);

	return vecY / length;
}

void Street::findIndicesLB(SPtr<Grid> grid, real initialSearchHeight)
{
	for (uint i = 0; i < numberOfCells; i++)
	{
		real x = getCoordinateX(i);
		real y = getCoordinateY(i);

		uint matrixIndex = grid->transCoordToIndex(x, y, initialSearchHeight);

		real xLB, yLB, zLB;
		grid->transIndexToCoords(matrixIndex, xLB, yLB, zLB);

		while (grid->getFieldEntry(matrixIndex) != BC_SOLID ||
			   grid->getFieldEntry(grid->transCoordToIndex(xLB, yLB, zLB-grid->getDelta())) != STOPPER_SOLID)
		{
			zLB -= grid->getDelta();
			matrixIndex = grid->transCoordToIndex(xLB, yLB, zLB);
		}

		std::stringstream msg;


		msg << "( " << x << ", " << y << " )" << "  ==>  ";
		msg << "( " << xLB << ", " << yLB << ", " << zLB << " ), type = [" << (int)grid->getFieldEntry(matrixIndex) << "], z = " << zLB << " \n";

		*logging::out << logging::Logger::INFO_LOW << msg.str();

		this->matrixIndicesLB.push_back(matrixIndex);
		this->sparseIndicesLB.push_back(grid->getSparseIndex(matrixIndex));
	}
}

void StreetPointFinder::prepareSimulationFileData()
{
	//////////////////////////////////////////////////////////////////////////
	// Concatenate sparseIndicesLB

	for (auto& street : this->streets) this->sparseIndicesLB.insert(this->sparseIndicesLB.end(), street.sparseIndicesLB.begin(), street.sparseIndicesLB.end());

	//////////////////////////////////////////////////////////////////////////
	// prepare vectors

	uint numberOfCells = this->sparseIndicesLB.size();

	mapNashToConc.resize(numberOfCells);

	std::vector<uint> indexMap(numberOfCells);
	std::iota(indexMap.begin(), indexMap.end(), 0);

	//////////////////////////////////////////////////////////////////////////
	// sort vectors

	std::stable_sort(indexMap.begin(),
		indexMap.end(),
		[&](uint lhs, uint rhs) {
		return this->sparseIndicesLB[lhs] <= this->sparseIndicesLB[rhs];
	});

	std::stable_sort(this->sparseIndicesLB.begin(),
		this->sparseIndicesLB.end(),
		[](uint lhs, uint rhs) {
		return lhs <= rhs;
	});
	//////////////////////////////////////////////////////////////////////////
	// invert idxMap

	{
		std::vector<uint> buffer = indexMap;
		for (uint idx = 0; idx < indexMap.size(); idx++)
			indexMap[buffer[idx]] = idx;
	}

	//////////////////////////////////////////////////////////////////////////
	// identify duplicates and find correct mapping indices

	std::vector<uint> reducedIndexMap(numberOfCells);

	uint currentSparseIndex = this->sparseIndicesLB[0];
	uint reducedIndex = 0;
	for (uint index = 1; index < numberOfCells; index++)
	{
		if (this->sparseIndicesLB[index] == currentSparseIndex)
		{
			reducedIndexMap[index] = reducedIndex;
		}
		else
		{
			currentSparseIndex = this->sparseIndicesLB[index];
			reducedIndexMap[index] = ++reducedIndex;
		}
	}

	for (uint index = 0; index < numberOfCells; index++)
	{
		mapNashToConc[index] = reducedIndexMap[indexMap[index]];
	}

	//////////////////////////////////////////////////////////////////////////
	// erase duplicated

	auto newEnd = std::unique(this->sparseIndicesLB.begin(), this->sparseIndicesLB.end());

	this->sparseIndicesLB.resize(std::distance(this->sparseIndicesLB.begin(), newEnd));

	//////////////////////////////////////////////////////////////////////////
}

void StreetPointFinder::readStreets(std::string filename)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readStreets( " << filename << " )" << "\n";

	uint numberOfStreets;

	std::ifstream file;

	file.open(filename.c_str());

	file >> numberOfStreets;

	for (uint i = 0; i < numberOfStreets; i++)
	{
		real xStart, yStart, xEnd, yEnd;

		real dx;

		file >> xStart >> yStart >> xEnd >> yEnd >> dx;

		streets.push_back(Street(xStart, yStart, xEnd, yEnd, dx));
	}

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::findIndicesLB(SPtr<Grid> grid, real initialSearchHeight)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::findIndicesLB()\n";

	for (auto& street : streets) street.findIndicesLB(grid, initialSearchHeight);

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeVTK(std::string filename, const std::vector<int>& cars)
{
	uint numberOfCells = 0;

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeVTK( " << filename << " )" << "\n";

	std::ofstream file;

	file.open(filename);

	prepareWriteVTK(file, numberOfCells);

	//////////////////////////////////////////////////////////////////////////

	file << "FIELD Label " << 3 << std::endl;

	//////////////////////////////////////////////////////////////////////////

	writeStreetsVTK(file, numberOfCells);

	writeLengthsVTK(file, numberOfCells);

	writeCarsVTK(file, numberOfCells, cars);

	////////////////////////////////////////////////////////////////////////////

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


void StreetPointFinder::writeReducedVTK(std::string filename, const std::vector<int>& cars)
{
	uint numberOfCells = 0;

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeVTK( " << filename << " )" << "\n";

	std::ofstream file;

	file.open(filename);

	prepareWriteVTK(file, numberOfCells);

	//////////////////////////////////////////////////////////////////////////

	file << "FIELD Label " << 1 << std::endl;

	//////////////////////////////////////////////////////////////////////////

	writeCarsVTK(file, numberOfCells, cars);

	////////////////////////////////////////////////////////////////////////////

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::prepareWriteVTK(std::ofstream & file, uint & numberOfCells)
{

	uint numberOfNodes = 0;

	for (auto& street : streets)
	{
		numberOfCells += street.numberOfCells;
		numberOfNodes += street.numberOfCells + 1;
	}	

	file << "# vtk DataFile Version 3.0\n";
	file << "by MeshGenerator\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	file << "POINTS " << numberOfNodes << " float" << std::endl;

	for (auto& street : streets)
	{
		for (uint i = 0; i <= street.numberOfCells; i++)
		{
			file << 0.5 * (street.getCoordinateX(i - 1) + street.getCoordinateX(i)) << " "
				<< 0.5 * (street.getCoordinateY(i - 1) + street.getCoordinateY(i)) << " " << 0.0 << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELLS " << numberOfCells << " " << 3 * numberOfCells << std::endl;


	uint nodeIndex = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			file << "2 " << nodeIndex << " " << nodeIndex + 1 << std::endl;
			nodeIndex++;
		}
		nodeIndex++;
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELL_TYPES " << numberOfCells << std::endl;

	for (uint i = 0; i < numberOfCells; i++) {
		file << "3" << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////

	file << "\nCELL_DATA " << numberOfCells << std::endl;
}


void StreetPointFinder::writeStreetsVTK(std::ofstream & file, uint numberOfCells)
{
	file << "StreetIndex 1 " << numberOfCells << " int" << std::endl;

	uint streetIndex = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			file << streetIndex << std::endl;
		}
		streetIndex++;
	}
}



void StreetPointFinder::writeCarsVTK(std::ofstream& file, uint numberOfCells, const std::vector<int>& cars)
{
	file << "Cars 1 " << numberOfCells << " float" << std::endl;

	uint index = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			if (index < cars.size())
				file << cars[index] << std::endl;
			else
				file << -1 << std::endl;
			index++;
		}
	}
}


void StreetPointFinder::writeLengthsVTK(std::ofstream & file, uint numberOfCells)
{
	file << "StreetLength 1 " << numberOfCells << " float" << std::endl;

	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			real length = std::sqrt((street.getCoordinateX(i) - street.getCoordinateX(0)) * (street.getCoordinateX(i) - street.getCoordinateX(0))
				+ (street.getCoordinateY(i) - street.getCoordinateY(0)) * (street.getCoordinateY(i) - street.getCoordinateY(0)));

			file << length << std::endl;
		}
	}
}


void StreetPointFinder::writeConnectionVTK(std::string filename, SPtr<Grid> grid)
{
	uint numberOfCells = 0;

	for (auto& street : streets)
	{
		numberOfCells += street.numberOfCells;
	}

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeConnectionVTK( " << filename << " )" << "\n";

	std::ofstream file;

	file.open(filename);

	file << "# vtk DataFile Version 3.0\n";
	file << "by MeshGenerator\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	file << "POINTS " << 2 * numberOfCells << " float" << std::endl;

	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			real xLB, yLB, zLB;
			grid->transIndexToCoords(street.matrixIndicesLB[i], xLB, yLB, zLB);

			file << street.getCoordinateX(i) << " " << street.getCoordinateY(i) << " " << 5.0 << std::endl;
			file << xLB << " " << yLB << " " << zLB << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELLS " << numberOfCells << " " << 3 * numberOfCells << std::endl;


	uint nodeIndex = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			file << "2 " << nodeIndex << " " << nodeIndex + 1 << std::endl;
			nodeIndex += 2;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELL_TYPES " << numberOfCells << std::endl;

	for (uint i = 0; i < numberOfCells; i++) {
		file << "3" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeSimulationFile(std::string gridPath, real concentration, uint numberOfLevels, uint level)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeSimulationFile( " << gridPath << "conc.dat )" << "\n";

	std::ofstream file;

	file.open(gridPath + "conc.dat");

	file << "concentration\n";

	file << numberOfLevels - 1 << "\n";

	for (uint currentLevel = 0; currentLevel < numberOfLevels; currentLevel++)
	{
		if (currentLevel == level)
		{
			uint numberOfCells = 0;
			for (auto& street : streets)
			{
				numberOfCells += street.numberOfCells;
			}

			file << numberOfCells << "\n";

			for (auto& street : streets)
			{
				for (auto& sparseIndexLB : street.sparseIndicesLB)
				{
					// + 1 for numbering shift between GridGenerator and VF_GPU
					file << sparseIndexLB + 1 << "\n";
				}
			}
		}
		else
		{
			file << "0\n";
		}
	}

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeStreetVectorFile(std::string gridPath, real concentration, uint numberOfLevels, uint level)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeStreetVectorFile( " << gridPath << "streetVector.dat )" << "\n";

	std::ofstream file;

	file.open(gridPath + "streetVector.dat");

	file << "streetVector\n";

	file << numberOfLevels - 1 << "\n";

	for (uint currentLevel = 0; currentLevel < numberOfLevels; currentLevel++)
	{
		if (currentLevel == level)
		{
			uint numberOfCells = 0;
			for (auto& street : streets)
			{
				numberOfCells += street.numberOfCells;
			}

			file << numberOfCells << "\n";

			for (auto& street : streets)
			{
				for (auto& sparseIndexLB : street.sparseIndicesLB)
				{
					// + 1 for numbering shift between GridGenerator and VF_GPU
					file << street.getVectorX() << " " << street.getVectorY() << "\n";
				}
			}
		}
		else
		{
			file << "0\n";
		}
	}

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeSimulationFileSorted(std::string gridPath, real concentration, uint numberOfLevels, uint level)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeSimulationFile( " << gridPath << "concSorted.dat )" << "\n";

	std::ofstream file;

	file.open(gridPath + "concSorted.dat");

	file << "concentration\n";

	file << numberOfLevels - 1 << "\n";

	for (uint currentLevel = 0; currentLevel < numberOfLevels; currentLevel++)
	{
		if (currentLevel == level)
		{
			file << this->sparseIndicesLB.size() << "\n";

			for (auto& sparseIndexLB : this->sparseIndicesLB)
			{
				// + 1 for numbering shift between GridGenerator and VF_GPU
				file << sparseIndexLB + 1 << "\n";
			}
		}
		else
		{
			file << "0\n";
		}
	}

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeMappingFile(std::string gridPath)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeMappingFile( " << gridPath << "mappingNashToConc.dat )" << "\n";

	std::ofstream file;

	file.open(gridPath + "mappingNashToConc.dat");

	file << this->mapNashToConc.size() << "\n";

	for (auto& index : this->mapNashToConc)
	{
		file << index << "\n";
	}

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Speed hackend by Stephan Lenz, not tested

void StreetPointFinder::write3DVTK(std::string filename, const std::vector<int>& cars)
{
	uint numberOfCells = 0;

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeVTK( " << filename << " )" << "\n";

	std::ofstream file;

	file.open(filename);

	prepareWrite3DVTK(file, numberOfCells, cars);

	////////////////////////////////////////////////////////////////////////////

	file.close();

	*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::prepareWrite3DVTK(std::ofstream & file, uint & numberOfCells, const std::vector<int>& cars)
{

	uint numberOfNodes = 0;

	for (auto& street : streets)
	{
		numberOfCells += street.numberOfCells;
		numberOfNodes += street.numberOfCells + 1;
	}

	file << "# vtk DataFile Version 3.0\n";
	file << "by MeshGenerator\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	uint index = 0;
	uint numberOfCars = 0;

	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			if (index < cars.size() && cars[index] != -1)
			{
				numberOfCars++;
			}

			index++;
		}
	}

	file << "POINTS " << 8 * numberOfCars << " float" << std::endl;

	index = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			if(index < cars.size() && cars[index] != -1 )
			{
				real xStart = 0.5 * (street.getCoordinateX(i - 1) + street.getCoordinateX(i));
				real yStart = 0.5 * (street.getCoordinateY(i - 1) + street.getCoordinateY(i));

				real xEnd = 0.5 * (street.getCoordinateX(i) + street.getCoordinateX(i + 1));
				real yEnd = 0.5 * (street.getCoordinateY(i) + street.getCoordinateY(i + 1));

				real vecX = xEnd - xStart;
				real vecY = yEnd - yStart;

				file << xStart + vecY << " " << yStart - vecX << " " << 0.0 << std::endl;
				file << xStart - vecY << " " << yStart + vecX << " " << 0.0 << std::endl;

				file << xEnd + vecY << " " << yEnd - vecX << " " << 0.0 << std::endl;
				file << xEnd - vecY << " " << yEnd + vecX << " " << 0.0 << std::endl;

				file << xStart + vecY << " " << yStart - vecX << " " << 1.5 << std::endl;
				file << xStart - vecY << " " << yStart + vecX << " " << 1.5 << std::endl;

				file << xEnd + vecY << " " << yEnd - vecX << " " << 1.5 << std::endl;
				file << xEnd - vecY << " " << yEnd + vecX << " " << 1.5 << std::endl;
			}

			index++;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELLS " << numberOfCars << " " << 9 * numberOfCars << std::endl;

	index = 0;
	uint carIndex = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			if (index < cars.size() && cars[index] != -1)
			{
				file << "8 " 
					 << 8 * carIndex + 0 << " "
					 << 8 * carIndex + 1 << " "
					 << 8 * carIndex + 3 << " "
					 << 8 * carIndex + 2 << " "
					 << 8 * carIndex + 4 << " "
					 << 8 * carIndex + 5 << " "
					 << 8 * carIndex + 7 << " "
					 << 8 * carIndex + 6 << " "
					 << std::endl;

				carIndex++;
			}
			index++;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	file << "CELL_TYPES " << numberOfCars << std::endl;

	for (uint i = 0; i < numberOfCars; i++) {
		file << "12" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////

	file << "\nCELL_DATA " << numberOfCars << std::endl;

	//////////////////////////////////////////////////////////////////////////

	file << "FIELD Label " << 1 << std::endl;

	//////////////////////////////////////////////////////////////////////////

	file << "Cars 1 " << numberOfCars << " float" << std::endl;

	index = 0;
	for (auto& street : streets)
	{
		for (uint i = 0; i < street.numberOfCells; i++)
		{
			if (index < cars.size() && cars[index] != -1)
				file << cars[index] << std::endl;
			
			index++;
		}
	}
}




