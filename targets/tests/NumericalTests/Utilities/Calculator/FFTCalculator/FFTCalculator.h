#ifndef FFTCALCULATOR_H
#define FFTCALCULATOR_H

#include "../Calculator.h"

#include <memory>
#include <vector>
#include <fftw3.h>

class Results;
class EvaluationParameter;
class TestResults;
class PhiAndNuTest;

class FFTCalculator : public Calculator
{
public:
	void calcAndCopyToTestResults();
	void setSimulationResults(std::shared_ptr<Results> simResults);

protected:
	FFTCalculator(double viscosity, std::shared_ptr<PhiAndNuTest> testResults);
	virtual void setVectorToCalc() = 0;

	std::shared_ptr<Results> simResults;
	std::vector<std::vector<double>> data;

private:
	void init();
	double calcNu();
	double calcNuDiff(double nu);
	double calcPhiDiff();
	std::vector<double> calcLinReg(std::vector<double> y);
	void calcLogAmplitudeForAllTimeSteps();
	void calcAmplitudeForAllTimeSteps();
	void calcPhiForAllTimeSteps();
	void calcFFT2D(unsigned int timeStep);
	void initDataForFFT(fftw_complex* input, unsigned int timeStep);
	void setFFTResults(fftw_complex* result, unsigned int timeStep);

	std::shared_ptr<PhiAndNuTest> testResults;
	std::vector<std::vector<double>> fftResultsIm;
	std::vector<std::vector<double>> fftResultsRe;
	std::vector<double> phi;
	std::vector<double> amplitude;
	std::vector<double> logAmplitude;

	bool fftCalculated;

	double lx, lz;
	double vis;
	double timeStepLength;
	int numberOfTimeSteps;

	double nu;
	double nudiff, phidiff;
};
#endif
