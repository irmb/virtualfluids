#include "D3Q27EsoTwist3DSplittedVectorEx.h"

D3Q27EsoTwist3DSplittedVectorEx::D3Q27EsoTwist3DSplittedVectorEx(int nx1, int nx2, int nx3, real value)
{
    this->NX1 = nx1;
    this->NX2 = nx2;
    this->NX3 = nx3;

    this->localDistributions = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
        new CbArray4D<real, IndexerX4X3X2X1>(13, nx1, nx2, nx3, value));
    this->nonLocalDistributions = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
        new CbArray4D<real, IndexerX4X3X2X1>(13, nx1, nx2, nx3, value));

    this->zeroDistributions =
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(nx1, nx2, nx3, value));
}
