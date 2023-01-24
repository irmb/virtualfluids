#include "RefineAroundGbObjectHelper.h"
#include <mpi/Communicator.h>
#include "OverlapBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>
#include <D3Q27TriFaceMeshInteractor.h>
#include <Grid3D.h>

RefineAroundGbObjectHelper::RefineAroundGbObjectHelper(SPtr<Grid3D> grid, int refineLevel,
                                                       SPtr<D3Q27TriFaceMeshInteractor> objectIter,
                                                       real startDistance, real stopDistance,
                                                       std::shared_ptr<vf::mpi::Communicator> comm)
    : grid(grid), refineLevel(refineLevel), objectIter(objectIter), startDistance(startDistance),
      stopDistance(stopDistance), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
RefineAroundGbObjectHelper::~RefineAroundGbObjectHelper(void) = default;
//////////////////////////////////////////////////////////////////////////
void RefineAroundGbObjectHelper::refine()
{
    using namespace vf::lbm::dir;

    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - start");

    int rank = grid->getRank();
    grid->setRank(0);

    objectIter->refineBlockGridToLevel(refineLevel, startDistance, stopDistance);

    RatioBlockVisitor ratioVisitor(refineLevel);
    grid->accept(ratioVisitor);

    RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
    grid->accept(ratioSmoothVisitor);

    OverlapBlockVisitor overlapVisitor(refineLevel, false);
    grid->accept(overlapVisitor);

    std::vector<int> dirs;
    for (int i = DIR_P00; i <= DIR_0MP; i++) {
        dirs.push_back(i);
    }
    SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
    grid->accept(interDirsVisitor);

    grid->setRank(rank);

    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - end");
}
