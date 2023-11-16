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
//! \file InitThixotropyBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef InitThixotropyBlockVisitor_H
#define InitThixotropyBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"

#include <muParser.h>

/*================================================================================*/
/*  D3Q27ETInitThixotropyBlockVisitor                                             */
/*                                                                                */




class InitThixotropyBlockVisitor : public Block3DVisitor
{
public:
    typedef std::numeric_limits<real> D3Q27RealLim;

public:
    InitThixotropyBlockVisitor();
    //D3Q27ETInitThixotropyBlockVisitor(LBMReal rho, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
    //! Constructor
    //! \param nu - viscosity
    //! \param rho - density
    //! \param vx1 - velocity in x
    //! \param vx2 - velocity in y
    //! \param vx3 - velocity in z
    //InitThixotropyBlockVisitor(LBMReal lambda /*LBMReal nu, LBMReal D, LBMReal rho, LBMReal vx1 = 0.0, LBMReal vx2 = 0.0, LBMReal vx3 = 0.0, LBMReal c=0.0, LBMReal f1 = 0.0, LBMReal f2 = 0.0, LBMReal f3 = 0.0*/);
    //////////////////////////////////////////////////////////////////////////
    //automatic vars are: x1,x2, x3
    //ussage example: setVx1("x1*0.01+x2*0.003")
    //////////////////////////////////////////////////////////////////////////
    //void setVx1(const mu::Parser& parser);
    //void setVx2(const mu::Parser& parser);
    //void setVx3(const mu::Parser& parser);
    //void setRho(const mu::Parser& parser);

    //void setVx1(const std::string& muParserString);
    //void setVx2(const std::string& muParserString);
    //void setVx3(const std::string& muParserString);
    //void setRho(const std::string& muParserString);
    ////////////////////////////////////////////////////////////////////////////
    //void setVx1(LBMReal vx1);
    //void setVx2(LBMReal vx2);
    //void setVx3(LBMReal vx3);
    //void setRho(LBMReal rho);
    //void setNu(LBMReal nu);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //void setf1(const mu::Parser& parser);
    //void setf2(const mu::Parser& parser);
    //void setf3(const mu::Parser& parser);
    void setLambda(const mu::Parser& parser);

    //void setf1(const std::string& muParserString);
    //void setf2(const std::string& muParserString);
    //void setf3(const std::string& muParserString);
    void setLambda(const std::string& muParserString);
    //////////////////////////////////////////////////////////////////////////
    //void setf1(LBMReal f1);
    //void setf2(LBMReal f2);
    //void setf3(LBMReal f3);
    void setLambda(real lambda);
    //void setD(LBMReal D);

    //void initialize(double* f, double x1, double x2, double x3, double vx1, double vx2, double vx3, double rho, UbTupleDouble3 coords, double dx, double o, bool NSE);

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

protected:
    void checkFunction(mu::Parser fct);
    typedef void(*CalcFeqsFct)(real* const& /*feq[27]*/, const real& /*(d)rho*/, const real& /*vx1*/, const real& /*vx2*/, const real& /*vx3*/);

private:
    mu::Parser muVx1;
    mu::Parser muVx2;
    mu::Parser muVx3;
    //mu::Parser muRho;
    //LBMReal nu;

    //mu::Parser muf1;
    //mu::Parser muf2;
    //mu::Parser muf3;
    mu::Parser muLambda;
    //LBMReal D;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H
