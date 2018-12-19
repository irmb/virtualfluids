#ifndef POST_PROCESSING_RESULTS_H
#define POST_PROCESSING_RESULTS_H

#include <vector>

class SimulationResults;

class PostProcessingResults
{
public:
	virtual void evaluate() = 0;

	virtual double getNuVx() = 0;
	virtual double getNuVy() = 0;
	virtual double getNuVz() = 0;
	virtual double getPhiDiffVx() = 0;
	virtual double getPhiDiffVy() = 0;
	virtual double getPhiDiffVz() = 0;

	virtual std::vector< double> getL2NormVx() = 0;
	virtual std::vector< double> getL2NormVy() = 0;
	virtual std::vector< double> getL2NormVz() = 0;
	virtual std::vector< double> getL2NormPress() = 0;
	virtual std::vector< double> getL2NormRho() = 0;

	virtual int getNumberOfXNodes() = 0;
	virtual int getNumberOfYNodes() = 0;
	virtual int getNumberOfZNodes() = 0;

private:


};
#endif