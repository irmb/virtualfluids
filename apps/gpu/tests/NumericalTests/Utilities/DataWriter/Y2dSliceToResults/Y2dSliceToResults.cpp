#include "Y2dSliceToResults.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "gpu/core/Parameter/Parameter.h"

#include <cmath>
#include <sstream>

std::shared_ptr<Y2dSliceToResults>
Y2dSliceToResults::getNewInstance(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo,
                                  unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults,
                                  unsigned int ySliceForCalculation)
{
    return std::shared_ptr<Y2dSliceToResults>(
        new Y2dSliceToResults(vectorWriterInfo, timeStepLength, simResults, ySliceForCalculation));
}

Y2dSliceToResults::Y2dSliceToResults()
{
}

Y2dSliceToResults::Y2dSliceToResults(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo,
                                     unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults,
                                     unsigned int ySliceForCalculation)
    : ToVectorWriter(vectorWriterInfo, timeStepLength)
{
    this->simResults = simResults;
    this->ySliceForCalculation = ySliceForCalculation;
}

void Y2dSliceToResults::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level)
{
    int timestep = t / timeStepLength;
    maxX = para->getGridX().at(level);
    maxY = para->getGridY().at(level);
    maxZ = para->getGridZ().at(level);

    int numberNodes = (maxX - 1) * (maxZ - 1);
    std::vector<double> x(numberNodes), y(numberNodes), z(numberNodes);
    std::vector<double> vx(numberNodes), vy(numberNodes), vz(numberNodes);
    std::vector<double> press(numberNodes), rho(numberNodes);
    std::vector<unsigned int> levels(numberNodes);

    ySliceForCalculation = maxY / 2;
    for (int posZ = 0; posZ < maxZ - 1; posZ++) {
        for (int posX = 0; posX < maxX - 1; posX++) {
            int posResults = CoordResults2DTo1D(posX, posZ);
            int posPara = CoordPara3DTo1D(posX, ySliceForCalculation, posZ);

            x.at(posResults) = (double)para->getParH(level)->coordinateX[posPara] - (double)1.0;
            y.at(posResults) = (double)para->getParH(level)->coordinateY[posPara] - (double)1.0;
            z.at(posResults) = (double)para->getParH(level)->coordinateZ[posPara] - (double)1.0;
            vx.at(posResults) = (double)para->getParH(level)->velocityX[posPara] * (double)para->getVelocityRatio();
            vy.at(posResults) = (double)para->getParH(level)->velocityY[posPara] * (double)para->getVelocityRatio();
            vz.at(posResults) = (double)para->getParH(level)->velocityZ[posPara] * (double)para->getVelocityRatio();
            press.at(posResults) = (double)para->getParH(level)->pressure[posPara] / (double)3.0 *
                                   (double)para->getDensityRatio() * (double)para->getVelocityRatio() *
                                   (double)para->getVelocityRatio();
            rho.at(posResults) = (double)para->getParH(level)->rho[posPara] / (double)3.0 *
                                 (double)para->getDensityRatio() * (double)para->getVelocityRatio() *
                                 (double)para->getVelocityRatio();
            levels.at(posResults) = level;
        }
    }
    simResults->addTimeStep(timestep, t, levels, x, y, z, vx, vy, vz, press, rho);
}

int Y2dSliceToResults::CoordPara3DTo1D(int x, int y, int z)
{
    return z * maxY * maxX + y * maxX + x + 1;
}

int Y2dSliceToResults::CoordResults2DTo1D(int x, int z)
{
    return z * (maxX - 1) + x;
}
