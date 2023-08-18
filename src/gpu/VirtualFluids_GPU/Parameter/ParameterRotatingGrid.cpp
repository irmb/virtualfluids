#include "ParameterRotatingGrid.h"
#include "Logger.h"
#include <iostream>

ParameterRotatingGrid::ParameterRotatingGrid(const std::array<real, 3> &centerPoint, const Axis &rotationalAxis,
                                             uint sizeOfNestedCoordinates)
    : rotationalAxis(rotationalAxis)
{
    parameterRotHost->centerPoint = centerPoint;
    parameterRotHost->sizeOfNestedCoordinates = sizeOfNestedCoordinates;
    parameterRotHost->memorySizeOfNestedCoordinates = sizeOfNestedCoordinates * sizeof(real);

    parameterRotDevice->centerPoint = centerPoint;
    parameterRotDevice->sizeOfNestedCoordinates = sizeOfNestedCoordinates;
    parameterRotDevice->memorySizeOfNestedCoordinates = sizeOfNestedCoordinates * sizeof(real);
}

void ParameterRotatingGrid::fillNestedCoordinateVectorsOnHost(const std::array<real *, 3> &globalCoordinates)
{
    if (parameterRotHost->nestedCoordinatesX == nullptr || parameterRotHost->nestedCoordinatesY == nullptr ||
        parameterRotHost->nestedCoordinatesZ == nullptr)
        throw std::runtime_error("Allocate host pointers of nestedCoordinates before filing them!");

    // #pragma omp parallel for
    for (uint index = 0; index < parameterRotHost->sizeOfNestedCoordinates; index++) {
        this->parameterRotHost->nestedCoordinatesX[index] = globalCoordinates[0][index] - parameterRotHost->centerPoint[0];
        this->parameterRotHost->nestedCoordinatesY[index] = globalCoordinates[1][index] - parameterRotHost->centerPoint[1];
        this->parameterRotHost->nestedCoordinatesZ[index] = globalCoordinates[2][index] - parameterRotHost->centerPoint[2];
    }
}

// void ParameterRotatingGrid::transformNestedToBase(const std::array<real *, 3> &globalCoordinates)
// {
// // #pragma omp parallel for
//     for (uint index = 0; index < nestedCoordinatesX.size(); index++) {
//         transformRotatingToGlobal(globalCoordinates[0][index], globalCoordinates[1][index], globalCoordinates[2][index],
//                                   nestedCoordinatesX[index], nestedCoordinatesY[index], nestedCoordinatesZ[index],
//                                   centerPoint[0], centerPoint[1], centerPoint[2], gridAngle[0], gridAngle[1], gridAngle[2]);
//     }
// }
