#include "SetUndefinedNodesBlockVisitor.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "ILBMKernel.h"

SetUndefinedNodesBlockVisitor::SetUndefinedNodesBlockVisitor(bool twoTypeOfConnectorsCheck)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), twoTypeOfConnectorsCheck(twoTypeOfConnectorsCheck)
{
}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    if (!block->hasInterpolationFlag())
        return;

    SPtr<ILBMKernel> kernel = block->getKernel();

    if (!kernel && (block->getRank() != grid->getRank()))
        return;

    // width of ghost layer
    // int gl = kernel->getGhostLayerWidth();
    int gl = 0;

    SPtr<BCArray3D> bcMatrix = kernel->getBCProcessor()->getBCArray();

    int minX1 = gl;
    int minX2 = gl;
    int minX3 = gl;

    int maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - gl;
    int maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - gl;
    int maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - gl;

    // int offset = 2;
    int offset = 3;

    if (block->hasInterpolationFlag(DIR_P00)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_M00)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0P0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0M0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_00P)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_00M)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PP0)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MM0)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PM0)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MP0)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_P0P)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_M0M)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_P0M)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_M0P)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0PP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0MM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0PM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_0MP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PPP)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MPP)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PMP)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MMP)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PPM)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MPM)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_PMM)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(DIR_MMM)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }

    //////////////////////////////////////////////////////////////////////////
    int offset2 = 1;
    int ll      = 0;

    minX1 = ll;
    minX2 = ll;
    minX3 = ll;

    maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - ll;
    maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - ll;
    maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - ll;

    if (block->hasInterpolationFlagFC(DIR_P00)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_M00)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0P0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0M0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_00P)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_00M)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PP0)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MM0)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PM0)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MP0)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_P0P)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_M0M)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_P0M)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_M0P)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0PP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0MM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0PM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_0MP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PPP)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MPP)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PMP)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MMP)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PPM)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MPM)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_PMM)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(DIR_MMM)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }

    // invert scaleCF blocks
    if (block->hasInterpolationFlagCF()) {
        if (block->hasInterpolationFlagFC() && twoTypeOfConnectorsCheck) {
            for (int i = DIR_P00; i <= DIR_MMM; i++) {
                UBLOG(logINFO, "FC in dir=" << i << " " << block->hasInterpolationFlagFC(i));
            }
            for (int i = DIR_P00; i <= DIR_MMM; i++) {
                UBLOG(logINFO, "CF in dir=" << i << " " << block->hasInterpolationFlagCF(i));
            }
            throw UbException(UB_EXARGS, "block " + block->toString() + " has CF and FC");
        }

        minX1 = gl;
        minX2 = gl;
        minX3 = gl;

        maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - gl;
        maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - gl;
        maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - gl;

        for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
            for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                    if (bcMatrix->isUndefined(ix1, ix2, ix3))
                        bcMatrix->setFluid(ix1, ix2, ix3);
                    else
                        bcMatrix->setUndefined(ix1, ix2, ix3);
                }
            }
        }

        return;
    }
}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::setNodesUndefined(int startix1, int endix1, int startix2, int endix2, int startix3,
                                                      int endix3, SPtr<BCArray3D> bcMatrix)
{
    for (int ix3 = startix3; ix3 <= endix3; ix3++)
        for (int ix2 = startix2; ix2 <= endix2; ix2++)
            for (int ix1 = startix1; ix1 <= endix1; ix1++) {
                bcMatrix->setUndefined(ix1, ix2, ix3);
            }
}
