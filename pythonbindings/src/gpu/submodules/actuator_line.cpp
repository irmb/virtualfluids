#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/ActuatorLine.h>
class PyActuatorLine : public ActuatorLine 
{
public:
    using ActuatorLine::ActuatorLine; // Inherit constructors
    void calcBladeForces() override 
    { 
        PYBIND11_OVERRIDE_NAME(void, ActuatorLine, "calc_blade_forces", calcBladeForces,); 
    }
};
namespace actuator_line
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        using arr = py::array_t<float, py::array::c_style>;
        
        py::class_<ActuatorLine, PreCollisionInteractor, PyActuatorLine, std::shared_ptr<ActuatorLine>>(parentModule, "ActuatorLine", py::dynamic_attr())
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
        .def_property("omega", &ActuatorLine::getOmega, &ActuatorLine::setOmega)
        .def_property("azimuth", &ActuatorLine::getAzimuth, &ActuatorLine::setAzimuth)
        .def_property("yaw", &ActuatorLine::getYaw, &ActuatorLine::setYaw)
        .def_property_readonly("n_blades", &ActuatorLine::getNBlades)
        .def_property_readonly("n_blade_nodes", &ActuatorLine::getNBladeNodes)
        .def_property_readonly("n_nodes", &ActuatorLine::getNNodes)
        .def_property_readonly("n_indices", &ActuatorLine::getNIndices)
        .def_property_readonly("density", &ActuatorLine::getDensity)
        .def_property_readonly("delta_t", &ActuatorLine::getDeltaT)
        .def_property_readonly("delta_x", &ActuatorLine::getDeltaX)
        .def_property_readonly("position_x", &ActuatorLine::getPositionX)
        .def_property_readonly("position_y", &ActuatorLine::getPositionY)
        .def_property_readonly("position_z", &ActuatorLine::getPositionZ)
        .def_property_readonly("position", [](ActuatorLine& al){ real position[3] = {al.getPositionX(), al.getPositionY(), al.getPositionZ()}; return arr(3, position); } )
        .def("get_radii", [](ActuatorLine& al){ return arr(al.getNBladeNodes(), al.getBladeRadii()); } )
        .def("get_blade_coords_x", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsX()); } )
        .def("get_blade_coords_y", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsY()); } )
        .def("get_blade_coords_z", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsZ()); } )        
        .def("get_blade_velocities_x", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesX()); } )
        .def("get_blade_velocities_y", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesY()); } )
        .def("get_blade_velocities_z", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesZ()); } )
        .def("get_blade_forces_x", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesX()); } )
        .def("get_blade_forces_y", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesY()); } )
        .def("get_blade_forces_z", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesZ()); } )
        .def("get_radii_device", [](ActuatorLine& al){ return arr(al.getNBladeNodes(), al.getBladeRadii()); } )
        .def("get_blade_coords_x_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsXD()); } )
        .def("get_blade_coords_y_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsYD()); } )
        .def("get_blade_coords_z_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeCoordsZD()); } )        
        .def("get_blade_velocities_x_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesXD()); } )
        .def("get_blade_velocities_y_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesYD()); } )
        .def("get_blade_velocities_z_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeVelocitiesZD()); } )
        .def("get_blade_forces_x_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesXD()); } )
        .def("get_blade_forces_y_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesYD()); } )
        .def("get_blade_forces_z_device", [](ActuatorLine& al){ return arr({al.getNBlades(), al.getNBladeNodes()}, al.getBladeForcesZD()); } )
        .def("set_preinit_blade_radii", [](ActuatorLine& al, arr radii){ al.setPreInitBladeRadii(static_cast<float *>(radii.request().ptr)); } )
        .def("set_blade_coords", [](ActuatorLine& al, arr coordsX, arr coordsY, arr coordsZ)
        { 
            al.setBladeCoords(static_cast<float *>(coordsX.request().ptr), static_cast<float *>(coordsY.request().ptr), static_cast<float *>(coordsZ.request().ptr)); 
        })
        .def("set_blade_velocities", [](ActuatorLine& al, arr velocitiesX, arr velocitiesY, arr velocitiesZ)
        { 
            al.setBladeVelocities(static_cast<float *>(velocitiesX.request().ptr), static_cast<float *>(velocitiesY.request().ptr), static_cast<float *>(velocitiesZ.request().ptr)); 
        })
        .def("set_blade_forces", [](ActuatorLine& al, arr forcesX, arr forcesY, arr forcesZ)
        { 
            al.setBladeForces(static_cast<float *>(forcesX.request().ptr), static_cast<float *>(forcesY.request().ptr), static_cast<float *>(forcesZ.request().ptr));
        })
        .def("set_blade_coords_device", [](ActuatorLine& al, arr coordsX, arr coordsY, arr coordsZ)
        { 
            al.setBladeCoordsD(static_cast<float *>(coordsX.request().ptr), static_cast<float *>(coordsY.request().ptr), static_cast<float *>(coordsZ.request().ptr)); 
        })
        .def("set_blade_velocities_device", [](ActuatorLine& al, arr velocitiesX, arr velocitiesY, arr velocitiesZ)
        { 
            al.setBladeVelocitiesD(static_cast<float *>(velocitiesX.request().ptr), static_cast<float *>(velocitiesY.request().ptr), static_cast<float *>(velocitiesZ.request().ptr)); 
        })
        .def("set_blade_forces_device", [](ActuatorLine& al, arr forcesX, arr forcesY, arr forcesZ)
        { 
            al.setBladeForcesD(static_cast<float *>(forcesX.request().ptr), static_cast<float *>(forcesY.request().ptr), static_cast<float *>(forcesZ.request().ptr)); 
        })
        .def("calc_blade_forces", &ActuatorLine::calcBladeForces);
    }
}