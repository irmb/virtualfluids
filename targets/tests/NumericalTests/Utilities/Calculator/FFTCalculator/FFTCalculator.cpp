#include "FFTCalculator.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Tests/PhiAndNuTest/PhiAndNuTest.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

void FFTCalculator::calc()
{
	init();

	nu = calcNu();
	nudiff = calcNuDiff(nu);
	phidiff = calcPhi();
}

double FFTCalculator::getNuDiff()
{
	return nudiff;
}

double FFTCalculator::getPhiDiff()
{
	return phidiff;
}

std::shared_ptr<FFTCalculator> FFTCalculator::getNewInstance(double viscosity)
{
	return std::shared_ptr<FFTCalculator>(new FFTCalculator(viscosity));
}

void FFTCalculator::setSimulationResults(std::shared_ptr<SimulationResults> simResults)
{
	this->simResults = simResults;
}

void FFTCalculator::setVectorToCalc(std::vector<std::vector<double>> data)
{
	this->data = data;
}

void FFTCalculator::init()
{
	fftResultsIm.clear();
	fftResultsRe.clear();
	phi.clear();
	amplitude.clear();
	logAmplitude.clear();
	lz = (double)simResults->getNumberOfZNodes();
	lx = (double)simResults->getNumberOfXNodes();
	timeStepLength = simResults->getTimeStepLength();
	numberOfTimeSteps = simResults->getNumberOfTimeSteps();
	fftCalculated = false;
}

FFTCalculator::FFTCalculator(double viscosity) : vis(viscosity)
{

}

double FFTCalculator::calcNu()
{
	calcLogAmplitudeForAllTimeSteps();
	std::vector<double> linReg = calcLinReg(logAmplitude);
	double nu = -(1.0 / (((2.0 * M_PI / lz) * (2.0 * M_PI / lz) + (2.0 * M_PI / lx)*(2.0 * M_PI / lx)) * timeStepLength)) * linReg.at(0);

	return nu;
}

double FFTCalculator::calcNuDiff(double nu)
{
	double nudiff = abs((nu - vis) / vis);
	return nudiff;
}

double FFTCalculator::calcPhi()
{
	phi = calcPhiForAllTimeSteps();
	std::vector<double> linReg = calcLinReg(phi);

	return linReg.at(0);
}

void FFTCalculator::calcLogAmplitudeForAllTimeSteps()
{
	calcAmplitudeForAllTimeSteps();
	logAmplitude.resize(amplitude.size());
	for(int tS=0; tS<amplitude.size();tS++)
		logAmplitude.at(tS) = log(amplitude.at(tS));
}

std::vector<double> FFTCalculator::calcLinReg(std::vector<double> y)
{
	std::vector<double> result;
	std::vector<double> x(y.size());
	double sumX = 0.0;
	double sumY = 0.0;

	for (int i = 0; i < y.size(); i++)
	{
		sumY += y.at(i);
		x.at(i) = i;
		sumX += i;
	}
	double avgX = sumX / y.size();
	double avgY = sumY / y.size();
	double zaehler = 0.0;
	double nenner = 0.0;
	for (int i = 0; i < y.size(); i++)
	{
		zaehler += (x.at(i) - avgX) * (y.at(i) - avgY);
		nenner += (x.at(i) - avgX) * (x.at(i) - avgX);
	}
	double a1 = zaehler / nenner;
	result.push_back(a1);
	double a0 = avgY - a1*avgX;
	result.push_back(a0);

	double ess = 0;
	double tss = 0;
	for (int i = 0; i < y.size(); i++)
	{
		ess += ((a0+a1*x.at(i))-avgY) * ((a0 + a1*x.at(i)) - avgY);
		tss += (y.at(i)-avgY) * (y.at(i) - avgY);
	}
	double r2 = ess / tss;
	result.push_back(r2);
	return result;
}

void FFTCalculator::calcAmplitudeForAllTimeSteps()
{
	if (fftCalculated == false) {
		for (int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
			calcFFT2D(timeStep);
		fftCalculated = true;
	}
	int pos = 2 + (lx - 1);
	for (int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
		amplitude.push_back(4.0 / (lx * lz)  * sqrt(fftResultsRe.at(timeStep).at(pos) * fftResultsRe.at(timeStep).at(pos) + fftResultsIm.at(timeStep).at(pos) * fftResultsIm.at(timeStep).at(pos)));
}

std::vector<double> FFTCalculator::calcPhiForAllTimeSteps()
{
	std::vector<double> result;

	if (fftCalculated == false) {
		for (int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
			calcFFT2D(timeStep);
		fftCalculated = true;
	}
	int pos = 2 + (lx - 1);
	for (int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
		result.push_back(atan(fftResultsIm.at(timeStep).at(pos) / fftResultsRe.at(timeStep).at(pos)));
	return result;
}

void FFTCalculator::calcFFT2D(unsigned int timeStep)
{
	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lx * lz);
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lx * lz);

	initDataForFFT(in, timeStep);

	fftw_plan p = fftw_plan_dft_2d(lz, lx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	setFFTResults(out, timeStep);

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}

void FFTCalculator::initDataForFFT(fftw_complex * input, unsigned int timeStep)
{
	for (int i = 0; i < data.at(timeStep).size(); i++)
	{
		input[i][0] = data.at(timeStep).at(i);
		input[i][1] = 0;
	}
}

void FFTCalculator::setFFTResults(fftw_complex * result, unsigned int timeStep)
{
	std::vector<double> fftRe, fftIm;
	fftRe.resize(data.at(timeStep).size());
	fftIm.resize(data.at(timeStep).size());

	for (int i = 0; i < data.at(timeStep).size(); i++)
	{
		fftRe.at(i) = result[i][0];
		fftIm.at(i) = result[i][1];
	}
	fftResultsIm.push_back(fftIm);
	fftResultsRe.push_back(fftRe);
}
