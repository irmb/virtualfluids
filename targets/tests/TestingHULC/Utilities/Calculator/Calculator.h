#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>
#include <vector>
#include <fftw3.h>

class Results;
class EvaluationParameter;

class Calulator {
public:
	Calulator(std::shared_ptr<Results> simResults, std::shared_ptr<EvaluationParameter> evaPara);


	double calcNu();
	double calcNuDiff(double nu);
	double calcPhiDiff();

private:
	std::vector<double> calcLinReg(std::vector<double> y);
	void calcLogAmplitudeForAllTimeSteps();
	void calcAmplitudeForAllTimeSteps();
	void calcPhiForAllTimeSteps();
	void calcFFT2D(unsigned int timeStep);
	void initDataForFFT(fftw_complex* input, unsigned int timeStep);
	void setFFTResults(fftw_complex* result, unsigned int timeStep);

	std::shared_ptr<Results> simResults;

	std::vector<std::vector<double>> data;
	std::vector<std::vector<double>> fftResultsIm;
	std::vector<std::vector<double>> fftResultsRe;
	std::vector<double> phi;
	std::vector<double> amplitude;
	std::vector<double> logAmplitude;

	bool fftCalculated;

	double lx, lz;
	double timeStepLength;
	double vis;
	int numberOfTimeSteps;
};
#endif
