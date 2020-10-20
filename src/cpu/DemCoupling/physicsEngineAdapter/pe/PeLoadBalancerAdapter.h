#ifndef PeLoadBalancerAdapter_h__
#define PeLoadBalancerAdapter_h__

#include "PointerDefinitions.h"
#include "blockforest/SetupBlockForest.h"

class Grid3D;
class Block3D;

class PeLoadBalancerAdapter
{
public:
    PeLoadBalancerAdapter(SPtr<Grid3D> grid, unsigned numberOfProcesses, int rank);
    walberla::uint_t operator()(walberla::SetupBlockForest &forest, const walberla::uint_t numberOfProcesses,
                                const walberla::memory_t perProcessMemoryLimit);
    unsigned getNumberOfProcesses() const { return numberOfProcesses; }
    int getRank() const { return rank; }

protected:
    SPtr<Block3D> getBlockByMinUniform(double minX1, double minX2, double minX3, SPtr<Grid3D> grid);

private:
    SPtr<Grid3D> grid;
    unsigned numberOfProcesses;
    int rank;
};

#endif // PeLoadBalancerAdapter_h__