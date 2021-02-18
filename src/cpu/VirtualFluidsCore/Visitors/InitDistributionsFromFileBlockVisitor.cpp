#include "InitDistributionsFromFileBlockVisitor.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "InitDensityLBMKernel.h"
#include "LBMKernel.h"
#include <basics/utilities/UbFileInputASCII.h>

InitDistributionsFromFileBlockVisitor::InitDistributionsFromFileBlockVisitor(/*LBMReal nu, */ LBMReal rho,
                                                                             std::string filename)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), /*nu(nu),*/ rho(rho)
{
    UbFileInputASCII in(filename);
    if (!in) {
        throw UbException(UB_EXARGS, "could not open file " + filename);
    }

    int nodesX1 = in.readInteger();
    int nodesX2 = in.readInteger();
    int nodesX3 = in.readInteger();

    matrix = CbArray4D<LBMReal, IndexerX4X3X2X1>(3, nodesX1, nodesX2, nodesX3, 0);

    for (int x3 = 0; x3 < nodesX3; x3++)
        for (int x2 = 0; x2 < nodesX2; x2++)
            for (int x1 = 0; x1 < nodesX1; x1++) {
                // for (int x1 = 0; x1<nodesX1; x1++)
                //   for (int x2 = 0; x2<nodesX2; x2++)
                //      for (int x3 = 0; x3<nodesX3; x3++)
                //      {
                matrix(Vx1, x1, x2, x3) = in.readDouble();
                matrix(Vx2, x1, x2, x3) = in.readDouble();
                matrix(Vx3, x1, x2, x3) = in.readDouble();
            }
}

//////////////////////////////////////////////////////////////////////////
void InitDistributionsFromFileBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace D3Q27System;

    if (!block)
        UB_THROW(UbException(UB_EXARGS, "block is not exist"));

    //   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    //   UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
    //   double dx = grid->getDeltaX(block);
    //   LBMReal o = LBMSystem::calcCollisionFactor(nu, block->getLevel());

    // Funktionszeiger
    typedef void (*CalcFeqsFct)(LBMReal *const & /*feq[27]*/, const LBMReal & /*(d)rho*/, const LBMReal & /*vx1*/,
                                const LBMReal & /*vx2*/, const LBMReal & /*vx3*/);
    CalcFeqsFct calcFeqsFct = NULL;

    LBMReal vx1, vx2, vx3;

    int gridRank  = grid->getRank();
    int blockRank = block->getRank();

    if (blockRank == gridRank && block->isActive()) {
        SPtr<ILBMKernel> kernel = block->getKernel();
        if (!kernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + block->toString());

        if (kernel->getCompressible())
            calcFeqsFct = &D3Q27System::calcCompFeq;
        else
            calcFeqsFct = &D3Q27System::calcIncompFeq;

        //      UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

        SPtr<BCArray3D> bcArray        = kernel->getBCProcessor()->getBCArray();
        SPtr<EsoTwist3D> distributions = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());

        LBMReal f[D3Q27System::ENDF + 1];

        //      size_t nx1 = distributions->getNX1();
        //      size_t nx2 = distributions->getNX2();
        //      size_t nx3 = distributions->getNX3();

        int minX1 = 0;
        int minX2 = 0;
        //      int minX3 = 0;

        int maxX1 = (int)bcArray->getNX1();
        int maxX2 = (int)bcArray->getNX2();
        int maxX3 = (int)bcArray->getNX3();

        int maxMX1 = (int)matrix.getNX2();
        int maxMX2 = (int)matrix.getNX3();
        int maxMX3 = (int)matrix.getNX4();

        int blockix1 = block->getX1();
        int blockix2 = block->getX2();
        int blockix3 = block->getX3();

        UbTupleInt3 blockNx = grid->getBlockNX();

        for (int ix3 = minX1; ix3 < maxX3; ix3++)
            for (int ix2 = minX2; ix2 < maxX2; ix2++)
                for (int ix1 = minX1; ix1 < maxX1; ix1++) {
                    int x1 = blockix1 * val<1>(blockNx) + ix1 - 1;
                    int x2 = blockix2 * val<2>(blockNx) + ix2 - 1;
                    int x3 = blockix3 * val<3>(blockNx) + ix3 - 1;

                    if (x1 == -1) {
                        x1 = maxMX1 - 1;
                    }
                    if (x2 == -1) {
                        x2 = maxMX2 - 1;
                    }
                    if (x3 == -1) {
                        x3 = maxMX3 - 1;
                    }

                    if (x1 == maxMX1) {
                        x1 = 1;
                    }
                    if (x2 == maxMX2) {
                        x2 = 1;
                    }
                    if (x3 == maxMX3) {
                        x3 = 1;
                    }

                    vx1 = matrix(Vx1, x1, x2, x3);
                    vx2 = matrix(Vx2, x1, x2, x3);
                    vx3 = matrix(Vx3, x1, x2, x3);

                    // int x1p, x2p, x3p;

                    ////x-derivative
                    // if (x1+1 >= maxMX1) x1p = x1;
                    // else  x1p = x1+1;
                    // double vx1Plusx1 = matrix(Vx1, x1p, x2, x3);
                    // double vx2Plusx1 = matrix(Vx2, x1p, x2, x3);
                    // double vx3Plusx1 = matrix(Vx3, x1p, x2, x3);

                    // if (x1-1 < minX1) x1p = x1;
                    // else  x1p = x1-1;
                    // double vx1Minusx1 = matrix(Vx1, x1p, x2, x3);
                    // double vx2Minusx1 = matrix(Vx2, x1p, x2, x3);
                    // double vx3Minusx1 = matrix(Vx3, x1p, x2, x3);

                    ////y-derivative
                    // if (x2+1 >= maxMX2) x2p = x2;
                    // else  x2p = x2+1;
                    // double vx1Plusx2 = matrix(Vx1, x1, x2p, x3);
                    // double vx2Plusx2 = matrix(Vx2, x1, x2p, x3);
                    // double vx3Plusx2 = matrix(Vx3, x1, x2p, x3);

                    // if (x2-1 < minX2) x2p = x2;
                    // else  x2p = x2-1;
                    // double vx1Minusx2 = matrix(Vx1, x1, x2p, x3);
                    // double vx2Minusx2 = matrix(Vx2, x1, x2p, x3);
                    // double vx3Minusx2 = matrix(Vx3, x1, x2p, x3);

                    ////z-derivative
                    // if (x3+1 >= maxMX3) x3p = x3;
                    // else  x3p = x3+1;
                    // double vx1Plusx3 = matrix(Vx1, x1, x2, x3p);
                    // double vx2Plusx3 = matrix(Vx2, x1, x2, x3p);
                    // double vx3Plusx3 = matrix(Vx3, x1, x2, x3p);

                    // if (x3-1 < minX3) x3p = x3;
                    // else  x3p = x3-1;
                    // double vx1Minusx3 = matrix(Vx1, x1, x2, x3);
                    // double vx2Minusx3 = matrix(Vx2, x1, x2, x3);
                    // double vx3Minusx3 = matrix(Vx3, x1, x2, x3);

                    // double ax = (vx1Plusx1 - vx1Minusx1) / (2.0);
                    // double bx = (vx2Plusx1 - vx2Minusx1) / (2.0);
                    // double cx = (vx3Plusx1 - vx3Minusx1) / (2.0);

                    // double ay = (vx1Plusx2 - vx1Minusx2) / (2.0);
                    // double by = (vx2Plusx2 - vx2Minusx2) / (2.0);
                    // double cy = (vx3Plusx2 - vx3Minusx2) / (2.0);

                    // double az = (vx1Plusx3 - vx1Minusx3) / (2.0);
                    // double bz = (vx2Plusx3 - vx2Minusx3) / (2.0);
                    // double cz = (vx3Plusx3 - vx3Minusx3) / (2.0);
                    // double eps_new = 1.0;
                    // LBMReal op = 1.;

                    // LBMReal feq[27];

                    // calcFeqsFct(feq, rho, vx1, vx2, vx3);

                    // double f_E = eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op) /
                    // (54.*o*op)); double f_N = f_E + eps_new *((2.*(ax - by)) / (9.*o)); double f_T = f_E + eps_new
                    // *((2.*(ax - cz)) / (9.*o)); double f_NE = eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op +
                    // ax*(5.*o + op) + by*(5.*o + op)) / (54.*o*op)); double f_SE = f_NE + eps_new *((ay + bx) /
                    // (9.*o)); double f_TE = eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op +
                    // ax*(5.*o + op)) / (54.*o*op)); double f_BE = f_TE + eps_new *((az + cx) / (9.*o)); double f_TN =
                    // eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op) /
                    // (54.*o*op)); double f_BN = f_TN + eps_new *((bz + cy) / (9.*o)); double f_ZERO = eps_new *((5.*(ax
                    // + by + cz)) / (9.*op)); double f_TNE = eps_new *(-(ay + az + bx + bz + cx + cy) / (72.*o)); double
                    // f_TSW = -eps_new *((ay + bx) / (36.*o)) - f_TNE; double f_TSE = -eps_new *((az + cx) / (36.*o)) -
                    // f_TNE; double f_TNW = -eps_new *((bz + cy) / (36.*o)) - f_TNE;

                    // f[E] = f_E + feq[E];
                    // f[W] = f_E + feq[W];
                    // f[N] = f_N + feq[N];
                    // f[S] = f_N + feq[S];
                    // f[T] = f_T + feq[T];
                    // f[B] = f_T + feq[B];
                    // f[NE] = f_NE + feq[NE];
                    // f[SW] = f_NE + feq[SW];
                    // f[SE] = f_SE + feq[SE];
                    // f[NW] = f_SE + feq[NW];
                    // f[TE] = f_TE + feq[TE];
                    // f[BW] = f_TE + feq[BW];
                    // f[BE] = f_BE + feq[BE];
                    // f[TW] = f_BE + feq[TW];
                    // f[TN] = f_TN + feq[TN];
                    // f[BS] = f_TN + feq[BS];
                    // f[BN] = f_BN + feq[BN];
                    // f[TS] = f_BN + feq[TS];
                    // f[TNE] = f_TNE + feq[TNE];
                    // f[TNW] = f_TNW + feq[TNW];
                    // f[TSE] = f_TSE + feq[TSE];
                    // f[TSW] = f_TSW + feq[TSW];
                    // f[BNE] = f_TSW + feq[BNE];
                    // f[BNW] = f_TSE + feq[BNW];
                    // f[BSE] = f_TNW + feq[BSE];
                    // f[BSW] = f_TNE + feq[BSW];
                    // f[REST] = f_ZERO + feq[REST];

                    calcFeqsFct(f, rho, vx1, vx2, vx3);

                    distributions->setDistribution(f, ix1, ix2, ix3);
                    distributions->setDistributionInv(f, ix1, ix2, ix3);
                    dynamicPointerCast<InitDensityLBMKernel>(kernel)->setVelocity(ix1, ix2, ix3, vx1, vx2, vx3);
                }
    }
}
//////////////////////////////////////////////////////////////////////////
