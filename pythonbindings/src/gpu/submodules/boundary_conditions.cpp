#include <pybind11/pybind11.h>
#include <gpu/GridGenerator/grid/BoundaryConditions/Side.h>
#include "gpu/VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"

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
        .def("set_velocity_boundary_condition", &BoundaryConditionFactory::setVelocityBoundaryCondition)
        .def("set_no_slip_boundary_condition", &BoundaryConditionFactory::setNoSlipBoundaryCondition)
        .def("set_slip_boundary_condition", &BoundaryConditionFactory::setSlipBoundaryCondition)
        .def("set_pressure_boundary_condition", &BoundaryConditionFactory::setPressureBoundaryCondition)
        .def("set_stress_boundary_condition", &BoundaryConditionFactory::setStressBoundaryCondition)
        .def("set_precursor_boundary_condition", &BoundaryConditionFactory::setPrecursorBoundaryCondition)
        .def("set_geometry_boundary_condition", &BoundaryConditionFactory::setGeometryBoundaryCondition);

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
        .value("PressureEquilibrium", BoundaryConditionFactory::PressureBC::PressureEquilibrium)
        .value("PressureEquilibrium2", BoundaryConditionFactory::PressureBC::PressureEquilibrium2)
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