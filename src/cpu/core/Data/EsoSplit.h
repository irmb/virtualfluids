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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Data Data
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef EsoSplit_h
#define EsoSplit_h

#include "D3Q27System.h"
#include "EsoTwist3D.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

//! \brief Class implements EsoTwist3D
//! \details EsoSplit uses three vectors to implement Esoteric Twist method
class EsoSplit : public EsoTwist3D
{
public:
    EsoSplit();
    //! \param nx1 number of nodes in x1 direction
    //! \param nx2 number of nodes in x2 direction
    //! \param nx3 number of nodes in x3 direction
    //! \param value initialisation value
    EsoSplit(size_t nx1, size_t nx2, size_t nx3, real value);
    //////////////////////////////////////////////////////////////////////////
    ~EsoSplit() override;
    //////////////////////////////////////////////////////////////////////////
    void swap() override;
    //////////////////////////////////////////////////////////////////////////
    void getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override;
    ////////////////////////////////////////////////////////////////////////
    void getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                     unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    real getPostCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    real getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX1() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX2() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX3() const override;
    //////////////////////////////////////////////////////////////////////////
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr getLocalDistributions();
    //////////////////////////////////////////////////////////////////////////
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr getNonLocalDistributions();
    //////////////////////////////////////////////////////////////////////////
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr getZeroDistributions();
    //////////////////////////////////////////////////////////////////////////
    void setNX1(size_t newNX1);
    void setNX2(size_t newNX2);
    void setNX3(size_t newNX3);
    void setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array);
    void setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array);
    void setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr array);

protected:
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;
    size_t NX1, NX2, NX3;

    friend class MPIIORestartSimulationObserver;
    friend class MPIIOMigrationSimulationObserver;
};

#endif

//! \}
