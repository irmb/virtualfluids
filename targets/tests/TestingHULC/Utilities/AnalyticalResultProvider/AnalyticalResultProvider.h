#ifndef ANALYTICAL_RESULT_H
#define ANALYTICAL_RESULT_H

#include <vector>
#include <memory>

class Results;
class Parameter;

class AnalyticalResultProvider
{
public:
	virtual std::vector < std::shared_ptr<Results> > getAnalyticalResults() = 0;

protected:
	virtual void calculate() = 0;
	virtual void init(std::vector< std::shared_ptr<Results> > simulationResults) = 0;

	std::vector< std::shared_ptr<Results> > analyticalResults;
	unsigned int numberNodes;
};
#endif