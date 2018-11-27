#include "FFTCalculator.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Tests/PhiAndNuTest/PhiAndNuTest.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

std::shared_ptr<FFTCalculator> FFTCalculator::getNewInstance(int lx, int lz, int timeStepLength)
{
	return std::shared_ptr<FFTCalculator>(new FFTCalculator(lx, lz, timeStepLength));
}

FFTCalculator::FFTCalculator(int lx, int lz, int timeStepLength)
{
	this->lx = (double)lx;
	this->lz = (double)lz;
	this->timeStepLength = (double)timeStepLength;
}

void FFTCalculator::calc(std::vector<std::vector<double>> data)
{
	this->data = data;
	init();

	nu = calcNu();
	phidiff = calcPhiDiff();
}

void FFTCalculator::init()
{
	fftResultsIm.clear();
	fftResultsRe.clear();
	fftCalculated = false;
}

double FFTCalculator::calcNu()
{
	std::vector<double> logAmplitude = calcLogAmplitudeForAllSteps();
	std::vector<double> linReg = calcLinReg(logAmplitude);
	double nu = -(1.0 / (((2.0 * M_PI / lz) * (2.0 * M_PI / lz) + (2.0 * M_PI / lx)*(2.0 * M_PI / lx)) * timeStepLength)) * linReg.at(0);

	return nu;
}

double FFTCalculator::calcPhiDiff()
{
	std::vector<double> phi = calcPhiForAllSteps();
	std::vector<double> linReg = calcLinReg(phi);

	return linReg.at(0);
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

std::vector<double> FFTCalculator::calcLogAmplitudeForAllSteps()
{
	std::vector<double> amplitude = calcAmplitudeForAllSteps();
	std::vector<double> logAmplitude;
	for (int i = 0; i < amplitude.size(); i++)
		logAmplitude.push_back(log(amplitude.at(i)));

	return logAmplitude;
}

std::vector<double> FFTCalculator::calcAmplitudeForAllSteps()
{
	std::vector<double> amplitude;
	if (fftCalculated == false) {
		for (int step = 0; step < data.size(); step++)
			calcFFT2D(step);
		fftCalculated = true;
	}
	int pos = 2 + (lx - 1);
	for (int step = 0; step < data.size(); step++)
		amplitude.push_back(4.0 / (lx * lz)  * sqrt(fftResultsRe.at(step).at(pos) * fftResultsRe.at(step).at(pos) + fftResultsIm.at(step).at(pos) * fftResultsIm.at(step).at(pos)));

	return amplitude;
}

std::vector<double> FFTCalculator::calcPhiForAllSteps()
{
	std::vector<double> phi;
	if (fftCalculated == false) {
		for (int step = 0; step < data.size(); step++)
			calcFFT2D(step);
		fftCalculated = true;
	}
	int pos = 2 + (lx - 1);
	for (int step = 0; step < data.size(); step++)
		phi.push_back(atan(fftResultsIm.at(step).at(pos) / fftResultsRe.at(step).at(pos)));
	return phi;
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

void FFTCalculator::initDataForFFT(fftw_complex * input, unsigned int step)
{

	for (int i = 0; i < data.at(step).size(); i++)
	{
		input[i][0] = data.at(step).at(i);
		input[i][1] = 0;
	}
}

void FFTCalculator::setFFTResults(fftw_complex * result, unsigned int step)
{
	std::vector<double> fftRe, fftIm;
	fftRe.resize(data.at(step).size());
	fftIm.resize(data.at(step).size());

	for (int i = 0; i < data.at(step).size(); i++)
	{
		fftRe.at(i) = result[i][0];
		fftIm.at(i) = result[i][1];
	}
	fftResultsIm.push_back(fftIm);
	fftResultsRe.push_back(fftRe);
}

double FFTCalculator::getNu()
{
	return nu;
}

double FFTCalculator::getPhiDiff()
{
	return phidiff;
}
