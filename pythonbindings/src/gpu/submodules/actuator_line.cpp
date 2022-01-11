#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/ActuatorLine.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace actuator_line
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<ActuatorLine, PreCollisionInteractor, std::shared_ptr<ActuatorLine>>(parentModule, "ActuatorLine")
        .def(py::init<  const uint,
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
                        "epsilon",
                        "turbine_pos_x", "turbine_pos_y", "turbine_pos_z", 
                        "diameter", 
                        "level", 
                        "delta_t", 
                        "delta_x")
        .def("get_blade_coords_x", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeCoordsX());})
        .def("get_blade_coords_y", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeCoordsY());})
        .def("get_blade_coords_z", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeCoordsZ());})        
        .def("get_blade_velocities_x", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeVelocitiesX());})
        .def("get_blade_velocities_y", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeVelocitiesY());})
        .def("get_blade_velocities_z", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeVelocitiesZ());})
        .def("get_blade_forces_x", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeForcesX());})
        .def("get_blade_forces_y", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeForcesY());})
        .def("get_blade_forces_z", [](ActuatorLine& al){ return py::array_t<float>(al.getNumberOfNodes(), al.getBladeForcesZ());})
        .def("set_blade_coords_x", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeCoordsX(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_coords_y", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeCoordsY(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_coords_z", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeCoordsZ(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_velocities_x", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeVelocitiesX(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_velocities_y", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeVelocitiesY(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_velocities_z", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeVelocitiesZ(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_forces_x", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeForcesX(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_forces_y", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeForcesY(static_cast<float *>(array.request().ptr)); })
        .def("set_blade_forces_z", [](ActuatorLine& al, py::array_t<float, py::array::c_style | py::array::forcecast> array){ al.setBladeForcesZ(static_cast<float *>(array.request().ptr)); })
        .def("calc_blade_forces", &ActuatorLine::calcBladeForces);
    }
}