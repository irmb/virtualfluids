#include "Y2dSliceToResults.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "../../Results/Results.h"
#include "Utilities/TestCondition/TestCondition.h"


Y2dSliceToResults::Y2dSliceToResults(std::shared_ptr<Results> simResults, unsigned int ySliceForCalculation, unsigned int startTimeY2dSliceToVector, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, std::shared_ptr<FileWriter> fileWriter, unsigned int startTimeDataWriter): ToVectorWriter(ySliceForCalculation, startTimeY2dSliceToVector, endTime, timeStepLength, writeFiles, fileWriter, startTimeDataWriter)
{
	this->simResults = simResults;
	counterTimeSteps = 0;
}

void Y2dSliceToResults::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level)
{
	counterTimeSteps++;

	maxX = para->getGridX().at(level);
	maxY = para->getGridY().at(level);
	maxZ = para->getGridZ().at(level);

	int numberNodes = (maxX - 1) * (maxZ - 1);
	std::vector<double> x(numberNodes), z(numberNodes);
	std::vector<double> vx(numberNodes), vz(numberNodes);
	std::vector<double> press(numberNodes), rho(numberNodes);

	for (int posZ = 0; posZ < maxZ - 1; posZ++)
	{
		for (int posX = 0; posX < maxX - 1; posX++)
		{
			int posResults = CoordResults2DTo1D(posX, posZ);
			int posPara = CoordPara3DTo1D(posX, ySliceForCalculation, posZ);

			x.at(posResults) = (double)para->getParH(level)->coordX_SP[posPara] - (double)1.0;
			z.at(posResults) = (double)para->getParH(level)->coordY_SP[posPara] - (double)1.0;
			vx.at(posResults) = (double)para->getParH(level)->vx_SP[posPara] * (double)para->getVelocityRatio();
			vz.at(posResults) = (double)para->getParH(level)->vz_SP[posPara] * (double)para->getVelocityRatio();
			press.at(posResults) = (double)para->getParH(level)->press_SP[posPara] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
			rho.at(posResults) = (double)para->getParH(level)->rho_SP[posPara] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
		}
	}
	simResults->addTimeStep(counterTimeSteps, t, x, z, vx, vz, press, rho);
	counterTimeSteps++;
}

int Y2dSliceToResults::CoordPara3DTo1D(int x, int y, int z)
{
	return z*maxY*maxX + y*maxX + x + 1;
}

int Y2dSliceToResults::CoordResults2DTo1D(int x, int z)
{
	return z * (maxX - 1) + x;
}
