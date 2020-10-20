//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file D3Q27EsoTwist3DSplittedVector.h
//! \ingroup Data
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef D3Q27EsoTwist3DSplittedVector_h
#define D3Q27EsoTwist3DSplittedVector_h

#include "D3Q27System.h"
#include "EsoTwist3D.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

//! \brief Class implements EsoTwist3D
//! \details D3Q27EsoTwist3DSplittedVector uses three vectors to implement Esoteric Twist method
class D3Q27EsoTwist3DSplittedVector : public EsoTwist3D
{
public:
    D3Q27EsoTwist3DSplittedVector();
    //! \param nx1 number of nodes in x1 direction
    //! \param nx2 number of nodes in x2 direction
    //! \param nx3 number of nodes in x3 direction
    //! \param value initialisation value
    D3Q27EsoTwist3DSplittedVector(size_t nx1, size_t nx2, size_t nx3, LBMReal value);
    //////////////////////////////////////////////////////////////////////////
    ~D3Q27EsoTwist3DSplittedVector() override;
    //////////////////////////////////////////////////////////////////////////
    void swap() override;
    //////////////////////////////////////////////////////////////////////////
    void getDistribution(LBMReal *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistribution(const LBMReal *const f, size_t x1, size_t x2, size_t x3) override;
    ////////////////////////////////////////////////////////////////////////
    void getDistributionInv(LBMReal *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistributionInv(const LBMReal *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistributionForDirection(const LBMReal *const f, size_t x1, size_t x2, size_t x3,
                                     unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistributionInvForDirection(const LBMReal *const f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX1() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX2() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX3() const override;
    //////////////////////////////////////////////////////////////////////////
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr getLocalDistributions();
    //////////////////////////////////////////////////////////////////////////
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr getNonLocalDistributions();
    //////////////////////////////////////////////////////////////////////////
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr getZeroDistributions();
    //////////////////////////////////////////////////////////////////////////
    void setNX1(size_t newNX1);
    void setNX2(size_t newNX2);
    void setNX3(size_t newNX3);
    void setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr array);
    void setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr array);
    void setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr array);

protected:
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;
    size_t NX1, NX2, NX3;

    friend class MPIIORestartCoProcessor;
    friend class MPIIOMigrationCoProcessor;
};

#endif
