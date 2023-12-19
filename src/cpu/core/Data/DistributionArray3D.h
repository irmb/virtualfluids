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

#ifndef DistributionArray3D_H
#define DistributionArray3D_H

#include <LBMSystem.h>

//! \brief Abstract class of data structure for LBM
class DistributionArray3D
{
public:
    DistributionArray3D() = default;

    virtual ~DistributionArray3D() = default;

    //! get number of nodes for x1 direction
    virtual size_t getNX1() const = 0;
    //! get number of nodes for x2 direction
    virtual size_t getNX2() const = 0;
    //! get number of nodes for x3 direction
    virtual size_t getNX3() const = 0;
    //! get distribution
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x2 coordinate x2
    //! \param x3 coordinate x3
    virtual void getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x2 coordinate x2
    //! \param x3 coordinate x3
    virtual void setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! get distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x2 coordinate x2
    //! \param x3 coordinate x3
    virtual void getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x1 coordinate x2
    //! \param x1 coordinate x3
    virtual void setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x1 coordinate x2
    //! \param x1 coordinate x3
    virtual void setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
    virtual void setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) = 0;
    virtual real getPostCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
    virtual void setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3, unsigned long int direction)  = 0;
    virtual void setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
    virtual real getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
    virtual void swap() = 0;

protected:
private:
};

#endif

//! \}
