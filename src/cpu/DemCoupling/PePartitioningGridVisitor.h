#ifndef PePartitioningGridVisitor_h
#define PePartitioningGridVisitor_h

#if defined VF_MPI

#include <PointerDefinitions.h>
#include <vector>

#include "Grid3DVisitor.h"

#include "PePhysicsEngineSolverAdapter.h"

#include <array>

////////////////////////////////////////////////////////////////////////
//! \brief The class implements domain decomposition with PE library
//! \author Konstantin Kutscher
//////////////////////////////////////////////////////////////////////////
namespace vf::mpi {class Communicator;}
class Grid3D;
class Block3D;
class DemCoProcessor;
// class walberla::blockforest::BlockForest;

class PePartitioningGridVisitor : public Grid3DVisitor
{
public:
    //! This describe different types of decomposition
    enum GraphType { LevelIntersected, LevelBased };

public:
    //! Constructor
    //! \param comm - communicator

    PePartitioningGridVisitor(std::shared_ptr<vf::mpi::Communicator> comm, std::shared_ptr<DemCoProcessor> dem);
    virtual ~PePartitioningGridVisitor();
    void visit(SPtr<Grid3D> grid) override;

protected:
    void collectData(SPtr<Grid3D> grid);
    void distributePartitionData(SPtr<Grid3D> grid);
    // void getBlocksByCuboid(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3,
    // std::vector<SPtr<Block3D>>& blocks, SPtr<Grid3D> grid);
    SPtr<Block3D> getBlockByMinUniform(double minX1, double minX2, double minX3, SPtr<Grid3D> grid);

private:
    std::shared_ptr<vf::mpi::Communicator> comm;
    std::shared_ptr<DemCoProcessor> dem;

    std::vector<int> ids;
    std::vector<int> ranks;

    std::shared_ptr<walberla::blockforest::BlockForest> forest;
};

#endif // VF_MPI
#endif
