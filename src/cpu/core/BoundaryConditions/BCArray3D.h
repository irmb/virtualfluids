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
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Sören Freudiger
//=======================================================================================

#ifndef BCArray_H
#define BCArray_H

#include "BoundaryConditions.h"
#include "CbArray3D.h"

#include <typeinfo>

#include <PointerDefinitions.h>

//! A class implements array to store boundary conditions flags
class BCArray3D
{
public:
    //////////////////////////////////////////////////////////////////////////
    BCArray3D();
    //////////////////////////////////////////////////////////////////////////
    BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3);
    //////////////////////////////////////////////////////////////////////////
    BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val);
    //////////////////////////////////////////////////////////////////////////
    virtual ~BCArray3D();
    //////////////////////////////////////////////////////////////////////////
    inline std::size_t getNX1() const;
    //////////////////////////////////////////////////////////////////////////
    inline std::size_t getNX2() const;
    //////////////////////////////////////////////////////////////////////////
    inline std::size_t getNX3() const;
    //////////////////////////////////////////////////////////////////////////
    void resize(std::size_t nx1, std::size_t nx2, std::size_t nx3);
    //////////////////////////////////////////////////////////////////////////
    void resize(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val);
    //////////////////////////////////////////////////////////////////////////
    bool validIndices(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline bool hasBC(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    void setBC(std::size_t x1, std::size_t x2, std::size_t x3, SPtr<BoundaryConditions> const &bc);
    //////////////////////////////////////////////////////////////////////////
    inline int getBCVectorIndex(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline const SPtr<BoundaryConditions> getBC(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline SPtr<BoundaryConditions> getBC(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    void setSolid(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    inline bool isSolid(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    void setFluid(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    // true : FLUID or BC
    // false: UNDEFINED or SOLID
    inline bool isFluid(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline bool isFluidWithoutBC(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline bool isUndefined(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    inline bool isUnvalidForCollision(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    void setUndefined(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    inline bool isInterfaceCF(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    void setInterfaceCF(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    inline bool isInterfaceFC(std::size_t x1, std::size_t x2, std::size_t x3) const;
    //////////////////////////////////////////////////////////////////////////
    void setInterfaceFC(std::size_t x1, std::size_t x2, std::size_t x3);
    //////////////////////////////////////////////////////////////////////////
    std::size_t getNumberOfSolidEntries() const;
    //////////////////////////////////////////////////////////////////////////
    std::size_t getNumberOfFluidEntries() const;
    //////////////////////////////////////////////////////////////////////////
    std::size_t getNumberOfFluidWithoutBCEntries() const;
    //////////////////////////////////////////////////////////////////////////
    std::size_t getNumberOfBCEntries() const;
    //////////////////////////////////////////////////////////////////////////
    std::size_t getNumberOfUndefinedEntries() const;
    //////////////////////////////////////////////////////////////////////////
    std::size_t getBCVectorSize() const;
    //////////////////////////////////////////////////////////////////////////
    std::string toString() const;
    //////////////////////////////////////////////////////////////////////////
    std::vector<int> &getBcindexmatrixDataVector();
    //////////////////////////////////////////////////////////////////////////
    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3, const int &ghostLayerWidth) const;

    static const int SOLID;
    static const int FLUID;
    static const int INTERFACECF;
    static const int INTERFACEFC;
    static const int UNDEFINED;

private:
    //////////////////////////////////////////////////////////////////////////
    void deleteBCAndSetType(std::size_t x1, std::size_t x2, std::size_t x3, int type);
    //////////////////////////////////////////////////////////////////////////
    void deleteBC(std::size_t x1, std::size_t x2, std::size_t x3);

    friend class MPIIORestartSimulationObserver;
    friend class MPIIOMigrationSimulationObserver;
    friend class MPIIOMigrationBESimulationObserver;

protected:
    //////////////////////////////////////////////////////////////////////////
    //-1 solid // -2 fluid -...
    CbArray3D<int, IndexerX3X2X1> bcindexmatrix;
    std::vector<SPtr<BoundaryConditions>> bcvector;
    std::vector<int> indexContainer;
};

//////////////////////////////////////////////////////////////////////////
inline std::size_t BCArray3D::getNX1() const { return bcindexmatrix.getNX1(); }
//////////////////////////////////////////////////////////////////////////
inline std::size_t BCArray3D::getNX2() const { return bcindexmatrix.getNX2(); }
//////////////////////////////////////////////////////////////////////////
inline std::size_t BCArray3D::getNX3() const { return bcindexmatrix.getNX3(); }
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::hasBC(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) >= 0;
}
//////////////////////////////////////////////////////////////////////////
inline int BCArray3D::getBCVectorIndex(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline const SPtr<BoundaryConditions> BCArray3D::getBC(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    int index = bcindexmatrix(x1, x2, x3);
    if (index < 0)
        return SPtr<BoundaryConditions>(); //=> NULL Pointer

    return bcvector[index];
}
//////////////////////////////////////////////////////////////////////////
inline SPtr<BoundaryConditions> BCArray3D::getBC(std::size_t x1, std::size_t x2, std::size_t x3)
{
    int index = bcindexmatrix(x1, x2, x3);
    if (index < 0)
        return SPtr<BoundaryConditions>(); //=> NULL Pointer

    return bcvector[index];
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isSolid(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) == SOLID;
}
//////////////////////////////////////////////////////////////////////////
// true : FLUID or BC
// false: UNDEFINED or SOLID
inline bool BCArray3D::isFluid(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    int tmp = bcindexmatrix(x1, x2, x3);
    return (tmp == FLUID || tmp >= 0);
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isFluidWithoutBC(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) == FLUID;
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isUndefined(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) == UNDEFINED;
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isUnvalidForCollision(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    const int type = bcindexmatrix(x1, x2, x3);
    return type == SOLID || type == UNDEFINED;
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isInterfaceCF(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) == INTERFACECF;
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isInterfaceFC(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    return bcindexmatrix(x1, x2, x3) == INTERFACEFC;
}
//////////////////////////////////////////////////////////////////////////
inline bool BCArray3D::isInsideOfDomain(const int &x1, const int &x2, const int &x3, const int &ghostLayerWidth) const
{
    const int minX1 = ghostLayerWidth;
    const int maxX1 = (int)this->getNX1() - 1 - ghostLayerWidth;
    const int minX2 = ghostLayerWidth;
    const int maxX2 = (int)this->getNX2() - 1 - ghostLayerWidth;
    const int minX3 = ghostLayerWidth;
    const int maxX3 = (int)this->getNX3() - 1 - ghostLayerWidth;

    return (!(x1 < minX1 || x1 > maxX1 || x2 < minX2 || x2 > maxX2 || x3 < minX3 || x3 > maxX3));
}

#endif

//! \}
