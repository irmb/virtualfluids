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
	static std::shared_ptr<FFTCalculator> getNewInstance(int lx, int lz, int timeStepLength);
	
	void calc(std::vector<std::vector<double>> data, bool transposeData);
	
	double getNu();
	double getPhiDiff();

private:
	FFTCalculator() {};
	FFTCalculator(int lx, int lz, int timeStepLength);
	void init();
	double calcNu();
	double calcPhiDiff();
	std::vector< double> calcPhiForAllSteps();
	std::vector< double> calcLinReg(std::vector<double> y);
	std::vector<double> calcLogAmplitudeForAllSteps();
	std::vector<double> calcAmplitudeForAllSteps();
	void calcFFT2D(unsigned int step);
	std::vector<std::vector<double>> transpose(std::vector<std::vector<double>>);
	void initDataForFFT(fftw_complex* input, unsigned int step);
	void setFFTResults(fftw_complex* result, unsigned int step);

	std::vector<std::vector<double>> data;
	std::vector<std::vector<double>> fftResultsIm;
	std::vector<std::vector<double>> fftResultsRe;
	
	bool fftCalculated;
	bool transposeData;
	double lx, lz;
	double timeStepLength;
	double nu;
	double phidiff;
};
#endif
