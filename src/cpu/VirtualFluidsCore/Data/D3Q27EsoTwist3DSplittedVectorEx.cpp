#include "D3Q27EsoTwist3DSplittedVectorEx.h"

D3Q27EsoTwist3DSplittedVectorEx::D3Q27EsoTwist3DSplittedVectorEx(int nx1, int nx2, int nx3, LBMReal value)
{
    this->NX1 = nx1;
    this->NX2 = nx2;
    this->NX3 = nx3;

    this->localDistributions = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
        new CbArray4D<LBMReal, IndexerX4X3X2X1>(13, nx1, nx2, nx3, value));
    this->nonLocalDistributions = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
        new CbArray4D<LBMReal, IndexerX4X3X2X1>(13, nx1, nx2, nx3, value));

    this->zeroDistributions =
        CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(nx1, nx2, nx3, value));
}
