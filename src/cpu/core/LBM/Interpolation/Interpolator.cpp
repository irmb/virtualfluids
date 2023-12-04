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
//! \file Interpolator.cpp
//! \ingroup Interpolation
//! \author Konstantin Kutscher
//=======================================================================================

#include "Interpolator.h"


struct Range
{
    Range(int maxX1, int maxX2, int maxX3) : m_maxX1(maxX1), m_maxX2(maxX2), m_maxX3(maxX3) {}
    inline bool operator()(int x1, int x2, int x3)
    {
        return x1 >= 0 && x1 < m_maxX1 && x2 >= 0 && x2 < m_maxX2 && x3 >= 0 && x3 < m_maxX3;
    }

    int m_maxX1;
    int m_maxX2;
    int m_maxX3;
};


//////////////////////////////////////////////////////////////////////////
void Interpolator::readICell(SPtr<DistributionArray3D> f, D3Q27ICell &icell, int x1, int x2, int x3)
{
    f->getPreCollisionDistribution(icell.BSW, x1, x2, x3);
    f->getPreCollisionDistribution(icell.BSE, x1 + 1, x2, x3);
    f->getPreCollisionDistribution(icell.BNW, x1, x2 + 1, x3);
    f->getPreCollisionDistribution(icell.BNE, x1 + 1, x2 + 1, x3);
    f->getPreCollisionDistribution(icell.TSW, x1, x2, x3 + 1);
    f->getPreCollisionDistribution(icell.TSE, x1 + 1, x2, x3 + 1);
    f->getPreCollisionDistribution(icell.TNW, x1, x2 + 1, x3 + 1);
    f->getPreCollisionDistribution(icell.TNE, x1 + 1, x2 + 1, x3 + 1);
}
//////////////////////////////////////////////////////////////////////////
void Interpolator::writeICell(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3)
{
    f->setPostCollisionDistribution(icell.BSW, x1, x2, x3);
    f->setPostCollisionDistribution(icell.BSE, x1 + 1, x2, x3);
    f->setPostCollisionDistribution(icell.BNW, x1, x2 + 1, x3);
    f->setPostCollisionDistribution(icell.BNE, x1 + 1, x2 + 1, x3);
    f->setPostCollisionDistribution(icell.TSW, x1, x2, x3 + 1);
    f->setPostCollisionDistribution(icell.TSE, x1 + 1, x2, x3 + 1);
    f->setPostCollisionDistribution(icell.TNW, x1, x2 + 1, x3 + 1);
    f->setPostCollisionDistribution(icell.TNE, x1 + 1, x2 + 1, x3 + 1);
}
//////////////////////////////////////////////////////////////////////////
void Interpolator::writeICellInv(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3)
{
    f->setPreCollisionDistribution(icell.BSW, x1, x2, x3);
    f->setPreCollisionDistribution(icell.BSE, x1 + 1, x2, x3);
    f->setPreCollisionDistribution(icell.BNW, x1, x2 + 1, x3);
    f->setPreCollisionDistribution(icell.BNE, x1 + 1, x2 + 1, x3);
    f->setPreCollisionDistribution(icell.TSW, x1, x2, x3 + 1);
    f->setPreCollisionDistribution(icell.TSE, x1 + 1, x2, x3 + 1);
    f->setPreCollisionDistribution(icell.TNW, x1, x2 + 1, x3 + 1);
    f->setPreCollisionDistribution(icell.TNE, x1 + 1, x2 + 1, x3 + 1);
}
//////////////////////////////////////////////////////////////////////////
void Interpolator::writeINode(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2, int x3)
{
    f->setPostCollisionDistribution(inode, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void Interpolator::writeINodeInv(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2,
                                           int x3)
{
    f->setPreCollisionDistribution(inode, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
bool Interpolator::iCellHasSolid(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3)
{
    for (int ix3 = x3; ix3 <= x3 + 1; ix3++)
        for (int ix2 = x2; ix2 <= x2 + 1; ix2++)
            for (int ix1 = x1; ix1 <= x1 + 1; ix1++) {
                if (bcArray->isSolid(ix1, ix2, ix3))
                    return true;
            }
    return false;
}
//////////////////////////////////////////////////////////////////////////
bool Interpolator::findNeighborICell(const SPtr<BCArray3D> bcArray, SPtr<DistributionArray3D> f,
                                               D3Q27ICell &icell, int maxX1, int maxX2, int maxX3, int x1, int x2,
                                               int x3, real &xoff, real &yoff, real &zoff)
{

    Range inRange(maxX1, maxX2, maxX3);

    // GoWest
    if (inRange(x1 - 1, x2, x3) && !iCellHasSolid(bcArray, x1 - 1, x2, x3)) {
        readICell(f, icell, x1 - 1, x2, x3);
        xoff = 1;
        yoff = 0;
        zoff = 0;
    }
    // GoEast
    else if (inRange(x1 + 2, x2, x3) &&
             !iCellHasSolid(bcArray, x1 + 1, x2, x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2, x3);
        xoff = -1;
        yoff = 0;
        zoff = 0;
    }
    // GoSouth
    else if (inRange(x1, x2 - 1, x3) && !iCellHasSolid(bcArray, x1, x2 - 1, x3)) {
        readICell(f, icell, x1, x2 - 1, x3);
        xoff = 0;
        yoff = 1;
        zoff = 0;
    }
    // GoNorth
    else if (inRange(x1, x2 + 2, x3) &&
             !iCellHasSolid(bcArray, x1, x2 + 1, x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1, x2 + 1, x3);
        xoff = 0;
        yoff = -1;
        zoff = 0;
    }
    // GoBottom
    else if (inRange(x1, x2, x3 - 1) && !iCellHasSolid(bcArray, x1, x2, x3 - 1)) {
        readICell(f, icell, x1, x2, x3 - 1);
        xoff = 0;
        yoff = 0;
        zoff = 1;
    }
    // GoTop
    else if (inRange(x1, x2, x3 + 2) &&
             !iCellHasSolid(bcArray, x1, x2, x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1, x2, x3 + 1);
        xoff = 0;
        yoff = 0;
        zoff = -1;
    }
    // GoNW
    else if (inRange(x1 - 1, x2 + 2, x3) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 + 1,
                            x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 + 1, x3);
        xoff = 1;
        yoff = -1;
        zoff = 0;
    }
    // GoNE
    else if (inRange(x1 + 2, x2 + 2, x3) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 + 1,
                            x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 + 1, x3);
        xoff = -1;
        yoff = -1;
        zoff = 0;
    }
    // GoSW
    else if (inRange(x1 - 1, x2 - 1, x3) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 - 1,
                            x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 - 1, x3);
        xoff = 1;
        yoff = 1;
        zoff = 0;
    }
    // GoSE
    else if (inRange(x1 + 2, x2 - 1, x3) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 - 1,
                            x3)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 - 1, x3);
        xoff = -1;
        yoff = 1;
        zoff = 0;
    }
    // GoBW
    else if (inRange(x1 - 1, x2, x3 - 1) && !iCellHasSolid(bcArray, x1 - 1, x2, x3 - 1)) {
        readICell(f, icell, x1 - 1, x2, x3 - 1);
        xoff = 1;
        yoff = 0;
        zoff = 1;
    }
    // GoBE
    else if (inRange(x1 + 2, x2, x3 - 1) &&
             !iCellHasSolid(bcArray, x1 + 1, x2,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2, x3 - 1);
        xoff = -1;
        yoff = 0;
        zoff = 1;
    }
    // GoBS
    else if (inRange(x1, x2 - 1, x3 - 1) && !iCellHasSolid(bcArray, x1, x2 - 1, x3 - 1)) {
        readICell(f, icell, x1, x2 - 1, x3 - 1);
        xoff = 0;
        yoff = 1;
        zoff = 1;
    }
    // GoBN
    else if (inRange(x1, x2 + 2, x3 - 1) &&
             !iCellHasSolid(bcArray, x1, x2 + 1,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1, x2 + 1, x3 - 1);
        xoff = 0;
        yoff = -1;
        zoff = 1;
    }
    // GoTW
    else if (inRange(x1 - 1, x2, x3 + 2) && !iCellHasSolid(bcArray, x1 - 1, x2, x3 + 1)) {
        readICell(f, icell, x1 - 1, x2, x3 + 1);
        xoff = 1;
        yoff = 0;
        zoff = -1;
    }
    // GoTE
    else if (inRange(x1 + 2, x2, x3 + 2) &&
             !iCellHasSolid(bcArray, x1 + 1, x2,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2, x3 + 1);
        xoff = -1;
        yoff = 0;
        zoff = -1;
    }
    // GoTS
    else if (inRange(x1, x2 - 1, x3 + 2) && !iCellHasSolid(bcArray, x1, x2 - 1, x3 + 1)) {
        readICell(f, icell, x1, x2 - 1, x3 + 1);
        xoff = 0;
        yoff = 1;
        zoff = -1;
    }
    // GoTN
    else if (inRange(x1, x2 + 2, x3 + 2) &&
             !iCellHasSolid(bcArray, x1, x2 + 1,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1, x2 + 1, x3 + 1);
        xoff = 0;
        yoff = -1;
        zoff = -1;
    }
    // GoTNW
    else if (inRange(x1 - 1, x2 + 2, x3 + 2) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 + 1,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 + 1, x3 + 1);
        xoff = 1;
        yoff = -1;
        zoff = -1;
    }
    // GoTNE
    else if (inRange(x1 + 2, x2 + 2, x3 + 2) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 + 1,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 + 1, x3 + 1);
        xoff = -1;
        yoff = -1;
        zoff = -1;
    }
    // GoTSE
    else if (inRange(x1 + 2, x2 - 1, x3 + 2) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 - 1,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 - 1, x3 + 1);
        xoff = -1;
        yoff = 1;
        zoff = -1;
    }
    // GoTSW
    else if (inRange(x1 - 1, x2 - 1, x3 + 2) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 - 1,
                            x3 + 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 - 1, x3 + 1);
        xoff = 1;
        yoff = 1;
        zoff = -1;
    }
    // GoBNW
    else if (inRange(x1 - 1, x2 + 2, x3 - 1) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 + 1,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 + 1, x3 - 1);
        xoff = 1;
        yoff = -1;
        zoff = 1;
    }
    // GoBNE
    else if (inRange(x1 + 2, x2 + 2, x3 - 1) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 + 1,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 + 1, x3 - 1);
        xoff = -1;
        yoff = -1;
        zoff = 1;
    }
    // GoBSE
    else if (inRange(x1 + 2, x2 - 1, x3 - 1) &&
             !iCellHasSolid(bcArray, x1 + 1, x2 - 1,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 + 1, x2 - 1, x3 - 1);
        xoff = -1;
        yoff = 1;
        zoff = 1;
    }
    // GoBSW
    else if (inRange(x1 - 1, x2 - 1, x3 - 1) &&
             !iCellHasSolid(bcArray, x1 - 1, x2 - 1,
                            x3 - 1)) // is next but one node in area (base node at 0,0,0)
    {
        readICell(f, icell, x1 - 1, x2 - 1, x3 - 1);
        xoff = 1;
        yoff = 1;
        zoff = 1;
    }
    // default
    else {
        // std::string err = "For x1="+StringUtil::toString(x1)+", x2=" + StringUtil::toString(x2)+", x3=" +
        // StringUtil::toString(x3)+
        //                  " interpolation is not implemented for other direction"+
        //                  " by using in: "+(std::string)typeid(*this).name()+
        //                  " or maybe you have a solid on the block boundary";
        // UB_THROW(UbException(UB_EXARGS, err));
        return false;
    }
    return true;
}
//////////////////////////////////////////////////////////////////////////
int Interpolator::iCellHowManySolids(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3)
{
    int count = 0;
    for (int ix3 = x3; ix3 <= x3 + 1; ix3++)
        for (int ix2 = x2; ix2 <= x2 + 1; ix2++)
            for (int ix1 = x1; ix1 <= x1 + 1; ix1++) {
                if (bcArray->isSolid(ix1, ix2, ix3))
                    count++;
            }
    return count;
}
