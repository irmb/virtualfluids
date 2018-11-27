#ifndef POST_PROCESSING_RESULTS_IMP_H
#define POST_PROCESSING_RESULTS_IMP_H

#include "PostProcessingResults.h"

#include <memory>
#include <string>

class AnalyticalResults;
class FFTCalculator;
class L2NormCalculator;
class SimulationResults;

class PostProcessingResultsImp : public PostProcessingResults
{
public:
	static std::shared_ptr< PostProcessingResultsImp> getNewInstance(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu, std::vector<std::string> dataToCalculateL2, unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu, unsigned int basicTimeStepL2Norm, unsigned int divergentTimeStepL2Norm);

	void evaluate();

	double getNuVx();
	double getNuVy();
	double getNuVz();
	double getPhiDiffVx();
	double getPhiDiffVy();
	double getPhiDiffVz();
	std::vector< double> getL2NormVx();
	std::vector< double> getL2NormVy();
	std::vector< double> getL2NormVz();
	std::vector< double> getL2NormPress();
	std::vector< double> getL2NormRho();
	int getNumberOfXNodes();
	int getNumberOfYNodes();
	int getNumberOfZNodes();

private:
	PostProcessingResultsImp() {};
	PostProcessingResultsImp(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu, std::vector<std::string> dataToCalculateL2, unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu, unsigned int basicTimeStepL2Norm, unsigned int divergentTimeStepL2Norm);

	std::vector<std::vector<double>> reduceDataToTimeSteps(std::vector<std::vector<double>> data, unsigned int startTimeStep, unsigned int endTimeStep);
	int calcTimeStepInResults(unsigned int timeStep);

	std::shared_ptr< SimulationResults> simResult;
	std::shared_ptr< FFTCalculator> fftCalculator;
	std::shared_ptr< L2NormCalculator> l2Normcalculator;
	std::shared_ptr< AnalyticalResults> analyticalResult;
	std::vector<std::string> dataToCalculatePhiAndNu;
	std::vector<std::string> dataToCalculateL2;
	unsigned int startTimeStepCalculationPhiNu;
	unsigned int endTimeStepCalculationPhiNu;
	unsigned int basicTimeStepL2Norm;
	unsigned int divergentTimeStepL2Norm;
	bool isEvaluated;

	double nuVx, nuVy, nuVz;
	double phiDiffVx, phiDiffVy, phiDiffVz;
	double l2VxBasic, l2VyBasic, l2VzBasic, l2RhoBasic, l2PressBasic;
	double l2VxDivergent, l2VyDivergent, l2VzDivergent, l2RhoDivergent, l2PressDivergent;
};
#endif