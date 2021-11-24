#include <pybind11/pybind11.h>
#include <VirtualFluids_GPU/PreCollisionInteractor/ActuatorLine.h>
#include <VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace actuator_line
{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module alModule = parentModule.def_submodule("actuator_line");

        py::class_<ActuatorLine, PreCollisionInteractor>(alModule, "ActuatorLine")
        .def(py::init<  const uint,
                        const uint,
                        const real,
                        const uint,
                        const real,
                        real, real, real,
                        const real,
                        int,
                        const real,
                        const real>(), 
                        "n_blades", 
                        "density", 
                        "n_blade_nodes", 
                        "turbine_pos_y", "turbine_pos_x", "turbine_pos_z", 
                        "diameter", 
                        "level", 
                        "delta_t", 
                        "delta_x")
        .def_property("blade_coordsX", &ActuatorLine::getBladeCoordsX, &ActuatorLine::setBladeCoordsX)
        .def_property("blade_coordsY", &ActuatorLine::getBladeCoordsY, &ActuatorLine::setBladeCoordsY)
        .def_property("blade_coordsZ", &ActuatorLine::getBladeCoordsZ, &ActuatorLine::setBladeCoordsZ)
        .def_property("blade_velocitiesX", &ActuatorLine::getBladeVelocitiesX, &ActuatorLine::setBladeVelocitiesX)
        .def_property("blade_velocitiesY", &ActuatorLine::getBladeVelocitiesY, &ActuatorLine::setBladeVelocitiesY)
        .def_property("blade_velocitiesZ", &ActuatorLine::getBladeVelocitiesZ, &ActuatorLine::setBladeVelocitiesZ)
        .def_property("blade_forcesX", &ActuatorLine::getBladeForcesX, &ActuatorLine::setBladeForcesX)
        .def_property("blade_forcesY", &ActuatorLine::getBladeForcesY, &ActuatorLine::setBladeForcesY)
        .def_property("blade_forcesZ", &ActuatorLine::getBladeForcesZ, &ActuatorLine::setBladeForcesZ)
;

        return alModule;
    }
    
    
}