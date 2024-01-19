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
//! \author Sven Marcus, Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbObject3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbLine3D.h>
#include <Interactors/Interactor3D.h>


namespace geometry
{
    namespace py = pybind11;

    template<class GeoObject>
    using py_geometry = py::class_<GeoObject, GbObject3D, std::shared_ptr<GeoObject>>;

    std::string GbPoint3D_repr_(const GbPoint3D &instance)
    {
        std::ostringstream stream;
        stream << "<GbPoint3D"
               << " x1: " << instance.getX1Coordinate()
               << " x2: " << instance.getX2Coordinate()
               << " x3: " << instance.getX3Coordinate() << ">";

        return stream.str();
    }

    void makeModule(py::module_ &parentModule)
    {
        py::module geometry = parentModule.def_submodule("geometry");

        py::class_<GbObject3D, std::shared_ptr<GbObject3D>>(geometry, "GbObject3D");

        py_geometry<GbPoint3D>(geometry, "GbPoint3D")
                .def(py::init())
                .def(py::init<double &, double &, double &>())
                .def(py::init<GbPoint3D *>())
                .def_property("x1", &GbPoint3D::getX1Coordinate, &GbPoint3D::setX1)
                .def_property("x2", &GbPoint3D::getX2Coordinate, &GbPoint3D::setX2)
                .def_property("x3", &GbPoint3D::getX3Coordinate, &GbPoint3D::setX3)
                .def("get_distance", &GbPoint3D::getDistance)
                .def("__repr__", &GbPoint3D_repr_);

        py_geometry<GbCuboid3D>(geometry, "GbCuboid3D")
                .def(py::init())
                .def(py::init<double &, double &, double &, double &, double &, double &>())
                .def(py::init<GbPoint3D *, GbPoint3D *>())
                .def(py::init<GbCuboid3D *>())
                .def_property("point1", &GbCuboid3D::getPoint1, &GbCuboid3D::setPoint1)
                .def_property("point2", &GbCuboid3D::getPoint2, &GbCuboid3D::setPoint2)
                .def("__repr__", [&](GbCuboid3D &instance)
                {
                    std::ostringstream stream;
                    stream << "<GbCuboid3D" << std::endl
                           << "point1: " << GbPoint3D_repr_(instance.getPoint1()) << std::endl
                           << "point2: " << GbPoint3D_repr_(instance.getPoint2()) << ">";
                    return stream.str();
                });

        py_geometry<GbLine3D>(geometry, "GbLine3D")
                .def(py::init())
                .def(py::init<GbPoint3D *, GbPoint3D *>())
                .def(py::init<GbLine3D>())
                .def_property("point1", &GbLine3D::getPoint1, &GbLine3D::setPoint1)
                .def_property("point2", &GbLine3D::getPoint2, &GbLine3D::setPoint2)
                .def("__repr__", [&](GbLine3D &instance)
                {
                    std::ostringstream stream;
                    stream << "<GbLine3D" << std::endl
                           << "point1: " << GbPoint3D_repr_(instance.getPoint1()) << std::endl
                           << "point2: " << GbPoint3D_repr_(instance.getPoint2()) << ">";
                    return stream.str();
                });


        py::class_<Interactor3D, std::shared_ptr<Interactor3D>>(geometry, "State")
                .def_readonly_static("SOLID", &Interactor3D::SOLID)
                .def_readonly_static("INVERSESOLID", &Interactor3D::INVERSESOLID)
                .def_readonly_static("TIMEDEPENDENT", &Interactor3D::TIMEDEPENDENT)
                .def_readonly_static("FLUID", &Interactor3D::FLUID)
                .def_readonly_static("MOVEABLE", &Interactor3D::MOVEABLE)
                .def_readonly_static("CHANGENOTNECESSARY", &Interactor3D::CHANGENOTNECESSARY);
    }

}