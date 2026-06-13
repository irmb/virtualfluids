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
//! \author Henry Korb
//=======================================================================================
#include "DataTypes.h"
#include "numpy_utils.h"
#include <basics/geometry3d/Axis.h>
#include <basics/geometry3d/GbSpatialData3D.h>
#include <memory>
#include <pybind11/detail/using_smart_holder.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <string>

namespace geometry3d
{
namespace py = pybind11;

template <typename T>
std::string getNameWithSuffix(const std::string& name);
template <>
std::string getNameWithSuffix<real>(const std::string& name)
{
    return name + "Real";
}
template <>
std::string getNameWithSuffix<real3>(const std::string& name)
{
    return name + "Real3";
}

template <typename T>
void makeSpatialDataSubmodules(py::module_& parentModule)
{
    using namespace numpy_utils;

    auto PyGbSpatialData3D =
        py::class_<GbSpatialData3D<T>, py::smart_holder>(parentModule, getNameWithSuffix<T>("GbSpatialData3D").c_str());

    py::native_enum<typename GbStructuredMesh3D<T>::InterpolationStrategy>(PyGbSpatialData3D, "InterpolationStrategy",
                                                                           "enum.Enum")
        .value("Trilinear", GbStructuredMesh3D<T>::InterpolationStrategy::Trilinear);

    py::class_<GbStructuredMesh3D<T>, GbSpatialData3D<T>, py::smart_holder>(
        parentModule, getNameWithSuffix<T>("GbStructuredMesh3D").c_str())
        .def(py::init([](py::array_t<real>& spacing, py::array_t<real>& origin, py::array_t<uint>& nPoints,
                         py::array_t<T>& values,
                         typename GbStructuredMesh3D<T>::InterpolationStrategy interpolationStrategy) {
                 return new GbStructuredMesh3D<T>(convertArrayToReal3(spacing.squeeze()),
                                                  convertArrayToReal3(origin.squeeze()),
                                                  convertArrayToArray<3, uint>(nPoints.squeeze()),
                                                  convertArrayToVec<T>(values.squeeze()), interpolationStrategy);
             }),
             py::arg("spacing"), py::arg("origin"), py::arg("n_points"), py::arg("values"),
             py::arg("interpolation_strategy"));

    auto pyGbPointCloud3D = py::class_<GbPointCloud3D<T>, GbSpatialData3D<T>, py::smart_holder>(
        parentModule, getNameWithSuffix<T>("GbPointCloud3D").c_str());

    (void)py::class_<typename GbPointCloud3D<T>::InterpolationStrategy, py::smart_holder>(pyGbPointCloud3D,
                                                                                          "InterpolationStrategy");

    py::class_<typename GbPointCloud3D<T>::InverseDistanceWeighing, typename GbPointCloud3D<T>::InterpolationStrategy,
               py::smart_holder>(pyGbPointCloud3D, "InverseDistanceWeighing")
        .def_static("make", GbPointCloud3D<T>::InverseDistanceWeighing::make, py::arg("n_points") = 10);

    py::class_<typename GbPointCloud3D<T>::NearestNeighbor, typename GbPointCloud3D<T>::InterpolationStrategy,
               py::smart_holder>(pyGbPointCloud3D, "NearestNeighbor")
        .def_static("make", GbPointCloud3D<T>::NearestNeighbor::make);

    pyGbPointCloud3D.def(
        py::init([](py::array_t<real>& coords, py::array_t<real>& data,
                    std::unique_ptr<typename GbPointCloud3D<T>::InterpolationStrategy> interpolationStrategy,
                    bool printTree) {
            return new GbPointCloud3D<T>(convertArrayToVecReal3(coords.squeeze()), convertArrayToVec<T>(data.squeeze()),
                                         std::move(interpolationStrategy), printTree);
        }),
        py::arg("points"), py::arg("data"), py::arg("interpolation_strategy"), py::arg("print_tree") = false);
}
void makeModule(py::module_& parentModule)
{
    using namespace numpy_utils;
    py::module geometry3dModule = parentModule.def_submodule("geometry3d");
    py::enum_<Axis>(geometry3dModule, "Axis").value("x", Axis::x).value("y", Axis::y).value("z", Axis::z);
    makeSpatialDataSubmodules<real>(geometry3dModule);
    makeSpatialDataSubmodules<real3>(geometry3dModule);
}
} // namespace geometry3d
