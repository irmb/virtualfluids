#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

#include "Sink/Sink.h"
#include "Source/Source.h"
#include "Junction/Junction.h"


class VF_PUBLIC CarDisplay {
public:
	CarDisplay(std::vector<int> **pcurrent, const uint safetyDistance);
	~CarDisplay() {};

	void initResults(uint timeSteps);

	void dispCurrentRoad() const;
	void dispResults(const std::vector<int> * neighbors, const std::vector<std::shared_ptr<Sink> > & sinks, const  std::vector<std::shared_ptr<Junction> > & junctions, const  std::vector<std::shared_ptr<Source> > & sources);
	void writeResultsToFile() const;

	void putCurrentIntoResults(uint step);

private:
	void visualizeSafetyDistanceForConsole(const std::vector<int> * neighbors);

	void dispJunctionsAtCell(uint index, const  std::vector<std::shared_ptr<Junction> > & junctions) const;
	void dispSinksAtCell(uint index, const std::vector<std::shared_ptr<Sink> > & sinks) const;
	void dispSourcesAtCell(uint index, const  std::vector<std::shared_ptr<Source> > & sources) const;

private:
	std::vector<std::vector<int> > results;		//saves the results of the calculation; x-axis = timesteps, y axis = positions
	
	std::vector<int> **ppcurrent;
	uint roadLength;
	const uint safetyDistance;

	uint timeSteps;
	uint writingStep;
};