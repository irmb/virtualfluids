#include "CarDisplay.h"

#include <fstream>
#include <iostream>
#include <iomanip>	//formatting output streams
#include <stdexcept>

#include "Utilities/VectorHelper.h"
#include "Utilities/safe_casting.h"
#include "Utilities/ConsoleColor.h"

CarDisplay::CarDisplay(std::vector<int> **pcurrent, const uint safetyDistance):
	safetyDistance{ safetyDistance }
{
	this->ppcurrent = pcurrent;
	roadLength = castSizeT_Uint((*pcurrent)->size());
}


void CarDisplay::initResults(uint timeSteps)
{
	this->timeSteps = timeSteps;

	results.resize(roadLength, std::vector<int>(1));

	for (uint i = 0; i < roadLength; i++) {
		results[i].resize(timeSteps + 1);
	}

	VectorHelper::fillVector(results, -1);
	putCurrentIntoResults(0);
}


void CarDisplay::putCurrentIntoResults(uint step)
{
	writingStep = step;
	for (uint i = 0; i < roadLength; i++) 
		results[i][writingStep] = (**ppcurrent)[i];	
}


void CarDisplay::writeResultsToFile() const
{
	try {


		std::fstream outFile("results.txt", std::fstream::out | std::fstream::trunc);
		if (outFile.is_open())
		{
			for (uint i = 0; i < results.size(); i++) {
				for (uint j = 0; j < results[i].size() - 1; j++)
					outFile << results[i][j] << " ";

				outFile << results[i][results[i].size() - 1];
				outFile << std::endl;
			}
			std::cout << "Finished writing data to file" << std::endl;
		}


		else
			throw std::runtime_error("Couldn't open file");

	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
	catch (...) {
		std::cerr << "unknown exception while writing to file" << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


void CarDisplay::dispCurrentRoad() const
{
	std::cout << "current: ( step: " << writingStep << " )" << std::endl;
	VectorHelper::dispVectorColour(**ppcurrent);
}


void CarDisplay::dispResults(const std::vector<int> * neighbors, const std::vector<std::shared_ptr<Sink> > & sinks, const  std::vector<std::shared_ptr<Junction> > & junctions, const  std::vector<std::shared_ptr<Source> > & sources)
{
	writeResultsToFile();

	visualizeSafetyDistanceForConsole(neighbors);
	reverse(results.begin(), results.end());

	for (uint i = 0; i < results.size(); i++) {

		dispJunctionsAtCell(i, junctions);
		dispSinksAtCell(i, sinks);
		dispSourcesAtCell(i, sources);

		for (uint j = 0; j < results[i].size(); j++) {
			VectorHelper::makeVectorOutputColourful(results[i][j]);
			std::cout << std::setw(4) << results[i][j];
		}

		std::cout << std::endl;
	}
	std::cout << std::endl;
	ConsoleColor::setDefaultWhite();
}


void CarDisplay::dispJunctionsAtCell(uint index, const  std::vector<std::shared_ptr<Junction> > & junctions)  const
{
	for (auto& junc : junctions) {
		ConsoleColor::setDefaultWhite();
		junc->dispJunction(index, roadLength);
	}
}


void CarDisplay::dispSinksAtCell(uint index, const std::vector<std::shared_ptr<Sink> > & sinks)  const
{
	for (auto& sink : sinks) {
		if (sink->getIndex() == roadLength - index - 1) {
			ConsoleColor::setBrightRed();
			std::cout << std::setw(4) << 1 - (sink->getPossibilityBeingBlocked());
			return;
		}
		std::cout << std::setw(4) << " ";
	}
}


void CarDisplay::dispSourcesAtCell(uint index, const  std::vector<std::shared_ptr<Source> > & sources)  const
{
	for (auto& source : sources) {
		if (source->getIndex() == roadLength - index - 1) {
			ConsoleColor::setBrightRed();
			std::cout << std::setw(4) << source->getPossibility();
			return;
		}
		std::cout << std::setw(4) << " ";
	}
}


void CarDisplay::visualizeSafetyDistanceForConsole(const std::vector<int>* neighbors)
{
	if (safetyDistance != 0) {
		int neighbor;
		for (uint step = 0; step <= timeSteps; step++) {
			for (uint i = 0; i < roadLength; i++) {
				if (results[i][step] > -1) {
					neighbor = (*neighbors)[i];
					for (uint j = 1; j <= safetyDistance; j++) {
						//junction or sink
						if (neighbor <= -1000)
							break;
						if (results[neighbor][step] > -1) {
							std::cerr << "safetyDistance was violated: timestep: " << step << "\t carIndex: " << i << std::endl;
							break;
						}
						else
							results[neighbor][step] = -5;
						neighbor = (*neighbors)[neighbor];
					}
				}
			}
		}
	}
}
