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

#ifndef ESOTWIST3D_H
#define ESOTWIST3D_H

#include "DistributionArray3D.h"
#include <LBMSystem.h>

//! \brief Abstract class for implementation of Esoteric Twist method
//! \details EsoTwist3D provide an interface for different implementations of Esoteric Twist method
//! <a href="https://doi.org/10.3390/computation5020019"><b>[ Geier et al., (2017), 10.3390/computation5020019]</b></a>
// Geier, M., & Schönherr, M. (2017). Esoteric twist: an efficient in-place streaming algorithmus for the lattice
// Boltzmann method on massively parallel hardware. Computation, 5(2), 19.

class EsoTwist3D : public DistributionArray3D
{
public:
    EsoTwist3D() = default;

    ~EsoTwist3D() override = default;

    //////////////////////////////////////////////////////////////////////////
    void swap() override = 0;
    //////////////////////////////////////////////////////////////////////////
    void getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override = 0;
    ////////////////////////////////////////////////////////////////////////
    void getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                     unsigned long int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    // virtual void getPostCollisionDistributionForDirection(real* const& f, const size_t& x1, const size_t& x2, const size_t&
    // x3, const unsigned long int& direction) = 0;
    //////////////////////////////////////////////////////////////////////////
    real getPostCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    real getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override = 0;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX1() const override = 0;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX2() const override = 0;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX3() const override = 0;
    //////////////////////////////////////////////////////////////////////////
};

#endif

//! \}
