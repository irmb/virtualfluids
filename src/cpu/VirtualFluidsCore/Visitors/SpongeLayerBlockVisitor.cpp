#include "SpongeLayerBlockVisitor.h"
#include "D3Q27System.h"
#include "LBMSystem.h"

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbException.h"
#include <geometry3d/GbCuboid3D.h>

using namespace std;

SpongeLayerBlockVisitor::SpongeLayerBlockVisitor(SPtr<GbCuboid3D> boundingBox, SPtr<LBMKernel> kernel, real nue,
                                                 int dir)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), boundingBox(boundingBox), kernel(kernel), nue(nue), dir(dir)
{
}
//////////////////////////////////////////////////////////////////////////
SpongeLayerBlockVisitor::~SpongeLayerBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    if (!boundingBox) {
        UB_THROW(UbException(UB_EXARGS, "The bounding box isn't set!"));
    }
    if (!kernel) {
        UB_THROW(UbException(UB_EXARGS, "The kernel isn't set!"));
    }
    if (kernel && (block->getRank() == grid->getRank())) {
        UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
        UbTupleDouble3 blockLengths = grid->getBlockLengths(block);

        real minX1 = val<1>(org);
        real minX2 = val<2>(org);
        real minX3 = val<3>(org);
        real maxX1 = val<1>(org) + val<1>(blockLengths);
        real maxX2 = val<2>(org) + val<2>(blockLengths);
        real maxX3 = val<3>(org) + val<3>(blockLengths);

        if (boundingBox->isCellInsideGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3)) {
            real collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
            kernel->setCollisionFactor(collFactor);
            kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
            kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
            kernel->setBlock(block);
            SPtr<LBMKernel> newKernel = kernel->clone();

            SPtr<DataSet3D> dataSet = block->getKernel()->getDataSet();
            if (!dataSet) {
                UB_THROW(UbException(UB_EXARGS,
                                     "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
            }

            newKernel->setDataSet(dataSet);

            SPtr<BCSet> bcProc = block->getKernel()->getBCSet();
            if (!bcProc) {
                UB_THROW(UbException(
                    UB_EXARGS, "It is not possible to change a BCSet in kernel! Old BCSet is not exist!"));
            }
            newKernel->setBCSet(bcProc);

            real oldCollFactor = newKernel->getCollisionFactor();

            UbTupleInt3 ixMin = grid->getBlockIndexes(boundingBox->getX1Minimum(), boundingBox->getX2Minimum(),
                                                      boundingBox->getX3Minimum());
            UbTupleInt3 ixMax = grid->getBlockIndexes(boundingBox->getX1Maximum(), boundingBox->getX2Maximum(),
                                                      boundingBox->getX3Maximum());

            real newCollFactor;

            if (dir == DIR_P00) {
                int ibX1      = block->getX1();
                int ibMax     = val<1>(ixMax) - val<1>(ixMin) + 1;
                real index  = (real)(ibX1 - val<1>(ixMin) + 1);
                newCollFactor = oldCollFactor - (oldCollFactor - c1o1) / (real)(ibMax)*index;
            } else if (dir == DIR_M00) {
                int ibX1      = block->getX1();
                int ibMax     = val<1>(ixMax) - val<1>(ixMin) + 1;
                real index  = (real)(ibX1 - val<1>(ixMin) + 1);
                newCollFactor = (oldCollFactor - c1o1) / (real)(ibMax)*index;
            } else if (dir == DIR_00P) {
                int ibX3      = block->getX3();
                int ibMax     = val<3>(ixMax) - val<3>(ixMin) + 1;
                real index  = (real)(ibX3 - val<3>(ixMin) + 1);
                newCollFactor = oldCollFactor - (oldCollFactor - c1o1) / (real)(ibMax)*index;
            } else if (dir == DIR_00M) {
                int ibX3      = block->getX3();
                int ibMax     = val<3>(ixMax) - val<3>(ixMin) + 1;
                real index  = (real)(ibX3 - val<3>(ixMin) + 1);
                newCollFactor = (oldCollFactor - c1o1) / (real)(ibMax)*index;
            } else
                UB_THROW(UbException(UB_EXARGS, "Problem: no orthogonal sponge layer!"));

            newKernel->setCollisionFactor(newCollFactor);
            block->setKernel(newKernel);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
