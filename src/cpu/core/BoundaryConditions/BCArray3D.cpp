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

#include "BCArray3D.h"

const int BCArray3D::SOLID       = -1;
const int BCArray3D::FLUID       = -2;
const int BCArray3D::INTERFACECF = -3;
const int BCArray3D::INTERFACEFC = -4;
const int BCArray3D::UNDEFINED   = -5;

//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D() = default;
//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3)
{
    bcindexmatrix.resize(nx1, nx2, nx3, UNDEFINED);
}
//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val)
{
    bcindexmatrix.resize(nx1, nx2, nx3, val);
}
//////////////////////////////////////////////////////////////////////////
BCArray3D::~BCArray3D() = default;
//////////////////////////////////////////////////////////////////////////
void BCArray3D::resize(std::size_t nx1, std::size_t nx2, std::size_t nx3) { bcindexmatrix.resize(nx1, nx2, nx3); }
//////////////////////////////////////////////////////////////////////////
void BCArray3D::resize(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val)
{
    bcindexmatrix.resize(nx1, nx2, nx3, val);
}
//////////////////////////////////////////////////////////////////////////
bool BCArray3D::validIndices(std::size_t x1, std::size_t x2, std::size_t x3) const
{
    if (x1 >= this->bcindexmatrix.getNX1())
        return false;
    if (x2 >= this->bcindexmatrix.getNX2())
        return false;
    if (x3 >= this->bcindexmatrix.getNX3())
        return false;
    return true;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setBC(std::size_t x1, std::size_t x2, std::size_t x3, SPtr<BoundaryConditions> const &bc)
{
    if (this->hasBC(x1, x2, x3)) {
        if (this->getBC(x1, x2, x3) == bc)
            return;
        else
            this->deleteBC(x1, x2, x3);
    }

    // if no vacant BCs available
    if (indexContainer.empty()) {
        bcvector.push_back(bc);
        bcindexmatrix(x1, x2, x3) = (int)bcvector.size() - 1;
    } else {
        int index                 = indexContainer.back();
        bcindexmatrix(x1, x2, x3) = index;
        bcvector[index]           = bc;
        indexContainer.pop_back();
    }
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setSolid(std::size_t x1, std::size_t x2, std::size_t x3)
{
    if (this->hasBC(x1, x2, x3))
        this->deleteBC(x1, x2, x3);
    bcindexmatrix(x1, x2, x3) = SOLID;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setFluid(std::size_t x1, std::size_t x2, std::size_t x3)
{
    if (this->hasBC(x1, x2, x3))
        this->deleteBC(x1, x2, x3);
    bcindexmatrix(x1, x2, x3) = FLUID;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setUndefined(std::size_t x1, std::size_t x2, std::size_t x3)
{
    if (this->hasBC(x1, x2, x3))
        this->deleteBC(x1, x2, x3);
    bcindexmatrix(x1, x2, x3) = UNDEFINED;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfSolidEntries() const
{
    const std::vector<int> &data = bcindexmatrix.getDataVector();
    std::size_t counter          = 0;
    for (std::size_t i = 0; i < data.size(); i++)
        if (data[i] == SOLID)
            counter++;
    return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfFluidEntries() const
{
    const std::vector<int> &data = bcindexmatrix.getDataVector();
    std::size_t counter          = 0;
    for (std::size_t i = 0; i < data.size(); i++) {
        int tmp = data[i];
        if (tmp == FLUID || tmp >= 0)
            counter++;
    }
    return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfFluidWithoutBCEntries() const
{
    const std::vector<int> &data = bcindexmatrix.getDataVector();
    std::size_t counter          = 0;
    for (std::size_t i = 0; i < data.size(); i++)
        if (data[i] == FLUID)
            counter++;
    return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfBCEntries() const
{
    const std::vector<int> &data = bcindexmatrix.getDataVector();
    std::size_t counter          = 0;
    for (std::size_t i = 0; i < data.size(); i++)
        if (data[i] >= 0)
            counter++;
    return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfUndefinedEntries() const
{
    const std::vector<int> &data = bcindexmatrix.getDataVector();
    std::size_t counter          = 0;
    for (std::size_t i = 0; i < data.size(); i++)
        if (data[i] == UNDEFINED)
            counter++;
    return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getBCVectorSize() const { return this->bcvector.size(); }
//////////////////////////////////////////////////////////////////////////
std::string BCArray3D::toString() const
{
    std::size_t solidCounter = 0;
    std::size_t fluidCounter = 0;
    std::size_t bcCounter    = 0;
    std::size_t undefCounter = 0;

    for (size_t x1 = 0; x1 < bcindexmatrix.getNX1(); x1++) {
        for (size_t x2 = 0; x2 < bcindexmatrix.getNX2(); x2++) {
            for (size_t x3 = 0; x3 < bcindexmatrix.getNX3(); x3++) {
                if (bcindexmatrix(x1, x2, x3) >= 0)
                    bcCounter++;
                else if (bcindexmatrix(x1, x2, x3) == FLUID)
                    fluidCounter++;
                else if (bcindexmatrix(x1, x2, x3) == SOLID)
                    solidCounter++;
                else if (bcindexmatrix(x1, x2, x3) == UNDEFINED)
                    undefCounter++;
                else
                    throw UbException(UB_EXARGS, "invalid matrixEntry");
            }
        }
    }

    std::size_t unrefEntriesInBcVector = 0;
    for (std::size_t i = 0; i < bcvector.size(); i++)
        if (!bcvector[i])
            unrefEntriesInBcVector++;

    std::stringstream text;
    text << "BCArray<" << typeid(SPtr<BoundaryConditions>).name() << "," << typeid(int).name() << ">";
    text << "[ entries: " << bcindexmatrix.getNX1() << "x" << bcindexmatrix.getNX2();
    text << "x" << bcindexmatrix.getNX3() << "=";
    text << bcindexmatrix.getNX1() * bcindexmatrix.getNX2() * bcindexmatrix.getNX3() << " ]:\n";
    text << " - #fluid entries : " << fluidCounter << std::endl;
    text << " - #bc    entries : " << bcCounter << std::endl;
    text << " - #solid entries : " << solidCounter << std::endl;
    text << " - #undef entries : " << undefCounter << std::endl;
    text << " - bcvector-entries      : " << bcvector.size() << " (empty ones: " << unrefEntriesInBcVector << ")\n";
    text << " - indexContainer-entries: " << indexContainer.size() << std::endl;

    return text.str();
}
//////////////////////////////////////////////////////////////////////////
std::vector<int> &BCArray3D::getBcindexmatrixDataVector() { return bcindexmatrix.getDataVector(); }

//////////////////////////////////////////////////////////////////////////
void BCArray3D::deleteBCAndSetType(std::size_t x1, std::size_t x2, std::size_t x3, int type)
{
    this->deleteBC(x1, x2, x3);

    // Assign matrix to new type
    bcindexmatrix(x1, x2, x3) = type;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::deleteBC(std::size_t x1, std::size_t x2, std::size_t x3)
{
    // check if BC exists at all
    int index = bcindexmatrix(x1, x2, x3);
    if (index < 0)
        return;

    // slide the released index into the index container
    indexContainer.push_back(index);

    //"delete" element
    bcvector[index] = SPtr<BoundaryConditions>();
}

//! \}
