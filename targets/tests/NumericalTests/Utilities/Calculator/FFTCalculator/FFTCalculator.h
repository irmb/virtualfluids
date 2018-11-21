#ifndef FFTCALCULATOR_H
#define FFTCALCULATOR_H

#include <memory>
#include <vector>
#include <fftw3.h>

class SimulationResults;
class EvaluationParameter;
class TestResults;
class PhiAndNuTest;

class FFTCalculator
{
public:
	static std::shared_ptr<FFTCalculator> getNewInstance(double viscosity);
	void setSimulationResults(std::shared_ptr<SimulationResults> simResults);
	void setVectorToCalc(std::vector<std::vector<double>> data);
	
	void calc();
	
	double getNuDiff();
	double getPhiDiff();

private:
	FFTCalculator(double viscosity);
	void init();
	double calcNu();
	double calcNuDiff(double nu);
	double calcPhi();
	std::vector< double> calcPhiForAllTimeSteps();
	std::vector< double> calcLinReg(std::vector<double> y);
	void calcLogAmplitudeForAllTimeSteps();
	void calcAmplitudeForAllTimeSteps();
	void calcFFT2D(unsigned int timeStep);
	void initDataForFFT(fftw_complex* input, unsigned int timeStep);
	void setFFTResults(fftw_complex* result, unsigned int timeStep);

	std::shared_ptr<SimulationResults> simResults;
	std::vector<std::vector<double>> data;
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
