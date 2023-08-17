#include "ParameterRotatingGrid.h"
#include "LBM/GPUHelperFunctions/CoordinateTransformation.h"

ParameterRotatingGrid::ParameterRotatingGrid(const std::array<real, 3> &centerPoint, const Axis &rotationalAxis)
    : centerPoint(centerPoint), rotationalAxis(rotationalAxis)
{
}

void ParameterRotatingGrid::initializeNestedCoordinates(const std::array<real *, 3> &globalCoordinates, uint numberOfNodes)
{
    this->nestedCoordinatesX.resize(numberOfNodes);
    this->nestedCoordinatesY.resize(numberOfNodes);
    this->nestedCoordinatesZ.resize(numberOfNodes);

// #pragma omp parallel for
    for (uint index = 0; index < numberOfNodes; index++) {
        this->nestedCoordinatesX[index] = globalCoordinates[0][index] - centerPoint[0];
        this->nestedCoordinatesY[index] = globalCoordinates[1][index] - centerPoint[1];
        this->nestedCoordinatesZ[index] = globalCoordinates[2][index] - centerPoint[2];
    }
}

void ParameterRotatingGrid::transformNestedToBase(const std::array<real *, 3> &globalCoordinates)
{
// #pragma omp parallel for
    for (uint index = 0; index < nestedCoordinatesX.size(); index++) {
        transformRotatingToGlobal(
            globalCoordinates[0][index], globalCoordinates[1][index], globalCoordinates[2][index], nestedCoordinatesX[index],
            nestedCoordinatesY[index], nestedCoordinatesZ[index], centerPoint[0], centerPoint[1], centerPoint[2],
            gridAngle * unitVectors.at(this->rotationalAxis)[0], gridAngle * unitVectors.at(this->rotationalAxis)[1],
            gridAngle * unitVectors.at(this->rotationalAxis)[2]);
    }
}
