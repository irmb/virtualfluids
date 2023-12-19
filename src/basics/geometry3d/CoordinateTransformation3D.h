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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef COORDINATETRANSFORMATION3D_H
#define COORDINATETRANSFORMATION3D_H

#include <cmath>
#include <sstream>
#include <string>

#include <basics/utilities/UbException.h>

#include <PointerDefinitions.h>

///////////////////////////////////////////////////////////////////////////////////////
//!
//! \brief A class provides 3d coordinate transformation
//! \details
//! description:     x1/x2/x3 = old, x1*/x2*/x3* = new
//!    x2
//!    ^             x*
//!    |            /
//!    |           2*
//!    4          /
//!    |         /
//!    3        1*                     => new coordsys is translated by originX1=originX2=originX3=2
//!    |       /                          new dx1=dx2=dx2=2 -> scaling by 2 in x1-,x2- und x3-direction
//!    2      /                           FIRST rotation by alpha around "x1" axis
//!    |       \                          THEN  rotation by beta  around "x2" axis
//!    1         \                        THEN  rotation by gamma around "x3" axis
//!    |           x1*
//!    |--1--2--3--4--5------------- > x1
//!
//!  Remark: It might be that the rotations around x1 and x3 axis are swapped.
//!
//////////////////////////////////////////////////////////////////////////////////////

class CoordinateTransformation3D
{
public:
    CoordinateTransformation3D();
    CoordinateTransformation3D(const double &originX1, const double &originX2, const double &originX3,
                               const double &dx1, const double &dx2, const double &dx3, const double &alpha,
                               const double &beta, const double &gamma);
    CoordinateTransformation3D(const double &originX1, const double &originX2, const double &originX3,
                               const double &dx1, const double &dx2, const double &dx3);
    CoordinateTransformation3D(CoordinateTransformation3D *transformation);

    void setTransformationValues(const double &originX1, const double &originX2, const double &originX3,
                                 const double &dx1, const double &dx2, const double &dx3, const double &alpha,
                                 const double &beta, const double &gamma);
    double getX1CoordinateOffset() const { return this->Tx1; } // Translation
    double getX2CoordinateOffset() const { return this->Tx2; }
    double getX3CoordinateOffset() const { return this->Tx3; }
    double getX1CoordinateScaling() const { return this->Sx1; } // Scaling
    double getX2CoordinateScaling() const { return this->Sx2; }
    double getX3CoordinateScaling() const { return this->Sx3; }
    double getRotationX1Angle() const { return this->alpha; }
    double getRotationX2Angle() const { return this->beta; }
    double getRotationX3Angle() const { return this->gamma; } // Rotation

    // Achtung die Winkel passen nicht ueberein -siehe setTransformationValues
    void setRotationX1Angle(double alpha)
    {
        this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, alpha,
                                      this->beta, this->gamma);
    }
    void setRotationX2Angle(double beta)
    {
        this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, this->alpha,
                                      beta, this->gamma);
    }
    void setRotationX3Angle(double gamma)
    {
        this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, this->alpha,
                                      this->beta, gamma);
    }

    void setActive(const bool &active);
    bool isActive() const { return this->active; }
    bool isTransformation() const { return this->transformation; }

    double transformForwardToX1Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformForwardToX2Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformForwardToX3Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformForwardToX1CoordinateIgnoringRotation(const double &x1) const;
    double transformForwardToX2CoordinateIgnoringRotation(const double &x2) const;
    double transformForwardToX3CoordinateIgnoringRotation(const double &x3) const;
    double transformBackwardToX1Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformBackwardToX2Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformBackwardToX3Coordinate(const double &x1, const double &x2, const double &x3) const;
    double transformBackwardToX1CoordinateIgnoringRotation(const double &x1) const;
    double transformBackwardToX2CoordinateIgnoringRotation(const double &x2) const;
    double transformBackwardToX3CoordinateIgnoringRotation(const double &x3) const;
    std::string toString() const;

private:
    double Tx1, Tx2, Tx3, Sx1, Sx2, Sx3, alpha, beta, gamma;

    double toX1factorX1, toX1factorX2, toX1factorX3, toX1delta;
    double toX2factorX1, toX2factorX2, toX2factorX3, toX2delta;
    double toX3factorX1, toX3factorX2, toX3factorX3, toX3delta;

    double fromX1factorX1, fromX1factorX2, fromX1factorX3, fromX1delta;
    double fromX2factorX1, fromX2factorX2, fromX2factorX3, fromX2delta;
    double fromX3factorX1, fromX3factorX2, fromX3factorX3, fromX3delta;

    bool active;
    bool transformation;

    friend class MPIIOSimulationObserver;
    friend class MPIIORestartSimulationObserver;
    friend class MPIIOMigrationSimulationObserver;
    friend class MPIIOMigrationBESimulationObserver;
    friend class CheckpointConverter;
};

#endif // COORDINATETRANSFORMATION3D_H

//! \}
