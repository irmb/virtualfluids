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
//! \file boindary_conditions.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <gpu/GridGenerator/grid/BoundaryConditions/Side.h>
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"

namespace boundary_conditions
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<SideType>(parentModule, "SideType")
        .value("MX", SideType::MX)
        .value("PX", SideType::PX)
        .value("MY", SideType::MY)
        .value("PY", SideType::PY)
        .value("MZ", SideType::MZ)
        .value("PZ", SideType::PZ)
        .value("GEOMETRY", SideType::GEOMETRY);

        py::class_<BoundaryConditionFactory>(parentModule, "BoundaryConditionFactory")
        .def(py::init<>())
        .def("set_velocity_boundary_condition", &BoundaryConditionFactory::setVelocityBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_no_slip_boundary_condition", &BoundaryConditionFactory::setNoSlipBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_slip_boundary_condition", &BoundaryConditionFactory::setSlipBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_pressure_boundary_condition", &BoundaryConditionFactory::setPressureBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_stress_boundary_condition", &BoundaryConditionFactory::setStressBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_precursor_boundary_condition", &BoundaryConditionFactory::setPrecursorBoundaryCondition, py::arg("boundary_condition_type"))
        .def("set_geometry_boundary_condition", &BoundaryConditionFactory::setGeometryBoundaryCondition, py::arg("boundary_condition_type"));

        py::enum_<BoundaryConditionFactory::VelocityBC>(parentModule, "VelocityBC")
        .value("VelocitySimpleBounceBackCompressible", BoundaryConditionFactory::VelocityBC::VelocitySimpleBounceBackCompressible)
        .value("VelocityIncompressible", BoundaryConditionFactory::VelocityBC::VelocityIncompressible)
        .value("VelocityCompressible", BoundaryConditionFactory::VelocityBC::VelocityCompressible)
        .value("VelocityAndPressureCompressible", BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible)
        .value("NotSpecified", BoundaryConditionFactory::VelocityBC::NotSpecified);


        py::enum_<BoundaryConditionFactory::NoSlipBC>(parentModule, "NoSlipBC")
        .value("NoSlipImplicitBounceBack", BoundaryConditionFactory::NoSlipBC::NoSlipImplicitBounceBack)
        .value("NoSlipBounceBack", BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack)
        .value("NoSlipIncompressible", BoundaryConditionFactory::NoSlipBC::NoSlipIncompressible)
        .value("NoSlipCompressible", BoundaryConditionFactory::NoSlipBC::NoSlipCompressible)
        .value("NoSlip3rdMomentsCompressible", BoundaryConditionFactory::NoSlipBC::NoSlip3rdMomentsCompressible);

        py::enum_<BoundaryConditionFactory::SlipBC>(parentModule, "SlipBC")
        .value("SlipIncompressible", BoundaryConditionFactory::SlipBC::SlipIncompressible)
        .value("SlipCompressible", BoundaryConditionFactory::SlipBC::SlipCompressible)
        .value("SlipBounceBack", BoundaryConditionFactory::SlipBC::SlipBounceBack)
        .value("SlipCompressibleTurbulentViscosity", BoundaryConditionFactory::SlipBC::SlipCompressibleTurbulentViscosity)
        .value("SlipPressureCompressibleTurbulentViscosity", BoundaryConditionFactory::SlipBC::SlipPressureCompressibleTurbulentViscosity)
        .value("NotSpecified", BoundaryConditionFactory::SlipBC::NotSpecified);

        py::enum_<BoundaryConditionFactory::PressureBC>(parentModule, "PressureBC")
        .value("PressureNonEquilibriumIncompressible", BoundaryConditionFactory::PressureBC::PressureNonEquilibriumIncompressible)
        .value("PressureNonEquilibriumCompressible", BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible)
        .value("OutflowNonReflective", BoundaryConditionFactory::PressureBC::OutflowNonReflective)
        .value("OutflowNonReflectivePressureCorrection", BoundaryConditionFactory::PressureBC::OutflowNonReflectivePressureCorrection)
        .value("NotSpecified", BoundaryConditionFactory::PressureBC::NotSpecified);

        py::enum_<BoundaryConditionFactory::StressBC>(parentModule, "StressBC")
        .value("StressCompressible", BoundaryConditionFactory::StressBC::StressCompressible)
        .value("StressBounceBack", BoundaryConditionFactory::StressBC::StressBounceBack)
        .value("StressPressureBounceBack", BoundaryConditionFactory::StressBC::StressPressureBounceBack)
        .value("NotSpecified", BoundaryConditionFactory::StressBC::NotSpecified);

        py::enum_<BoundaryConditionFactory::PrecursorBC>(parentModule, "PrecursorBC")
        .value("VelocityPrecursor", BoundaryConditionFactory::PrecursorBC::VelocityPrecursor)
        .value("DistributionsPrecursor", BoundaryConditionFactory::PrecursorBC::DistributionsPrecursor)
        .value("NotSpecified", BoundaryConditionFactory::PrecursorBC::NotSpecified);
    }
}