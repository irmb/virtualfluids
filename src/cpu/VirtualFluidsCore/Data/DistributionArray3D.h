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
//! \file DistributionArray3D.h
//! \ingroup Data
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
    virtual void getDistribution(real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x2 coordinate x2
    //! \param x3 coordinate x3
    virtual void setDistribution(const real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! get distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x2 coordinate x2
    //! \param x3 coordinate x3
    virtual void getDistributionInv(real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x1 coordinate x2
    //! \param x1 coordinate x3
    virtual void setDistributionInv(const real *const f, size_t x1, size_t x2, size_t x3) = 0;
    //! set distribution in inverse order
    //! \param f distribution
    //! \param x1 coordinate x1
    //! \param x1 coordinate x2
    //! \param x1 coordinate x3
    virtual void setDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                             unsigned long int direction)                               = 0;
    virtual void setDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) = 0;
    virtual real getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction)      = 0;
    virtual void setDistributionInvForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                unsigned long int direction)                            = 0;
    virtual void setDistributionInvForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                unsigned long int direction)                            = 0;
    virtual real getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction)         = 0;
    virtual void swap()                                                                                 = 0;

protected:
private:
};

#endif
