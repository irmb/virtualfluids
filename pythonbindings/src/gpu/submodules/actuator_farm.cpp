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
//! \file actuator_farm.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/ActuatorFarm.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>
class PyActuatorFarm : public ActuatorFarm 
{
public:
    using ActuatorFarm::ActuatorFarm; // Inherit constructors
    void calcBladeForces() override 
    { 
        PYBIND11_OVERRIDE_NAME(void, ActuatorFarm, "calc_blade_forces", calcBladeForces); 
    }
};
namespace actuator_farm
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        using arr = py::array_t<float, py::array::c_style>;
        
        py::class_<ActuatorFarm, PreCollisionInteractor, PyActuatorFarm, std::shared_ptr<ActuatorFarm>>(parentModule, "ActuatorFarm", py::dynamic_attr())
        .def(py::init<  const uint,
                        const real,
                        const uint,
                        const real,
                        int,
                        const real,
                        const real,
                        const bool>(), 
                        py::arg("number_of_blades_per_turbine"), 
                        py::arg("density"), 
                        py::arg("number_of_nodes_per_blade"), 
                        py::arg("epsilon"),
                        py::arg("level"), 
                        py::arg("delta_t"), 
                        py::arg("delta_x"),
                        py::arg("use_host_arrays"))
        .def_property_readonly("number_of_turbines", &ActuatorFarm::getNumberOfTurbines)
        .def_property_readonly("number_of_nodes_per_blade", &ActuatorFarm::getNumberOfNodesPerBlade)
        .def_property_readonly("number_of_blades_per_turbine", &ActuatorFarm::getNumberOfBladesPerTurbine)
        .def_property_readonly("number_of_grid_nodes", &ActuatorFarm::getNumberOfGridNodes)
        .def_property_readonly("number_of_indices", &ActuatorFarm::getNumberOfIndices)
        .def_property_readonly("density", &ActuatorFarm::getDensity)
        .def_property_readonly("delta_t", &ActuatorFarm::getDeltaT)
        .def_property_readonly("delta_x", &ActuatorFarm::getDeltaX)

        .def("add_turbine", &ActuatorFarm::addTurbine, py::arg("posX"), py::arg("posY"), py::arg("posZ"), py::arg("diameter"), py::arg("omega"), py::arg("azimuth"), py::arg("yaw"), py::arg("bladeRadii"))

        .def("get_turbine_pos", [](ActuatorFarm& al, uint turbine){ real position[3] = {al.getTurbinePosX(turbine), al.getTurbinePosY(turbine), al.getTurbinePosZ(turbine)}; return arr(3,  position); }, py::arg("turbine"))
        .def("get_turbine_azimuth", &ActuatorFarm::getTurbineAzimuth, py::arg("turbine"))
        .def("get_turbine_yaw", &ActuatorFarm::getTurbineYaw, py::arg("turbine"))
        .def("get_turbine_omega", &ActuatorFarm::getTurbineOmega, py::arg("turbine"))
        .def("get_all_azimuths", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllAzimuths()); } )
        .def("get_all_yaws", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllYaws()); } )
        .def("get_all_omegas", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllOmegas()); } )
        .def("get_all_turbine_pos_x", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosX()); } )
        .def("get_all_turbine_pos_y", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosY()); } )
        .def("get_all_turbine_pos_z", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosZ()); } )
    
        .def("get_all_blade_radii", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfNodesPerBlade()}, al.getAllBladeRadii()); } )
        .def("get_all_blade_coords_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsX()); } )
        .def("get_all_blade_coords_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsY()); } )
        .def("get_all_blade_coords_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsZ()); } )        
        .def("get_all_blade_velocities_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesX()); } )
        .def("get_all_blade_velocities_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesY()); } )
        .def("get_all_blade_velocities_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesZ()); } )
        .def("get_all_blade_forces_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesX()); } )
        .def("get_all_blade_forces_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesY()); } )
        .def("get_all_blade_forces_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesZ()); } )

        .def("get_turbine_blade_radii", [](ActuatorFarm& al, uint turbine){ return arr(al.getNumberOfNodesPerBlade(), al.getTurbineBladeRadiiDevice(turbine)); } , py::arg("turbine"))
        .def("get_turbine_blade_coords_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsZDevice(turbine)); }, py::arg("turbine") )        
        .def("get_turbine_blade_velocities_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesZDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesZDevice(turbine)); }, py::arg("turbine") )

        .def("get_all_blade_radii_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t>(al.getAllBladeRadiiDevice()); } )
        .def("get_all_blade_coords_x_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeCoordsXDevice()); } )
        .def("get_all_blade_coords_y_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeCoordsYDevice()); } )
        .def("get_all_blade_coords_z_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeCoordsZDevice()); } )        
        .def("get_all_blade_velocities_x_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeVelocitiesXDevice()); } )
        .def("get_all_blade_velocities_y_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeVelocitiesYDevice()); } )
        .def("get_all_blade_velocities_z_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeVelocitiesZDevice()); } )
        .def("get_all_blade_forces_x_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeForcesXDevice()); } )
        .def("get_all_blade_forces_y_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeForcesYDevice()); } )
        .def("get_all_blade_forces_z_device", [](ActuatorFarm& al) -> intptr_t { return reinterpret_cast<intptr_t> (al.getAllBladeForcesZDevice()); } )

        .def("get_turbine_blade_radii_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeRadiiDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeCoordsXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeCoordsYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeCoordsZDevice(turbine)); }, py::arg("turbine") )        
        .def("get_turbine_blade_velocities_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeVelocitiesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeVelocitiesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeVelocitiesZDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeForcesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeForcesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return reinterpret_cast<intptr_t>(al.getTurbineBladeForcesZDevice(turbine)); }, py::arg("turbine") )

        .def("set_all_azimuths", [](ActuatorFarm& al, arr azimuths){ al.setAllAzimuths(static_cast<float *>(azimuths.request().ptr)); }, py::arg("azimuths"))
        .def("set_all_yaws", [](ActuatorFarm& al, arr yaws){ al.setAllYaws(static_cast<float *>(yaws.request().ptr)); }, py::arg("yaws"))
        .def("set_all_omegas", [](ActuatorFarm& al, arr omegas){ al.setAllOmegas(static_cast<float *>(omegas.request().ptr)); }, py::arg("omegas"))

        .def("set_turbine_azimuth", &ActuatorFarm::setTurbineAzimuth, py::arg("turbine"), py::arg("azimuth"))
        .def("set_turbine_yaw", &ActuatorFarm::setTurbineYaw, py::arg("turbine"), py::arg("yaw"))
        .def("set_turbine_omega", &ActuatorFarm::setTurbineOmega, py::arg("turbine"), py::arg("omega"))

        .def("set_all_blade_coords", [](ActuatorFarm& al, arr coordsX, arr coordsY, arr coordsZ)
        { 
            al.setAllBladeCoords(static_cast<float *>(coordsX.request().ptr), static_cast<float *>(coordsY.request().ptr), static_cast<float *>(coordsZ.request().ptr)); 
        }, py::arg("blade_coords_x"), py::arg("blade_coords_y"), py::arg("blade_coords_z") )
        .def("set_all_blade_velocities", [](ActuatorFarm& al, arr velocitiesX, arr velocitiesY, arr velocitiesZ)
        { 
            al.setAllBladeVelocities(static_cast<float *>(velocitiesX.request().ptr), static_cast<float *>(velocitiesY.request().ptr), static_cast<float *>(velocitiesZ.request().ptr)); 
        }, py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z") )
        .def("set_all_blade_forces", [](ActuatorFarm& al, arr forcesX, arr forcesY, arr forcesZ)
        { 
            al.setAllBladeForces(static_cast<float *>(forcesX.request().ptr), static_cast<float *>(forcesY.request().ptr), static_cast<float *>(forcesZ.request().ptr));
        }, py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z") )     
        .def("set_turbine_blade_coords", [](ActuatorFarm& al, uint turbine, arr coordsX, arr coordsY, arr coordsZ)
        { 
            al.setTurbineBladeCoords(turbine, static_cast<float *>(coordsX.request().ptr), static_cast<float *>(coordsY.request().ptr), static_cast<float *>(coordsZ.request().ptr)); 
        }, py::arg("turbine"), py::arg("blade_coords_x"), py::arg("blade_coords_y"), py::arg("blade_coords_z") )
        .def("set_turbine_blade_velocities", [](ActuatorFarm& al, uint turbine, arr velocitiesX, arr velocitiesY, arr velocitiesZ)
        {
            al.setTurbineBladeVelocities(turbine, static_cast<float *>(velocitiesX.request().ptr), static_cast<float *>(velocitiesY.request().ptr), static_cast<float *>(velocitiesZ.request().ptr)); 
        }, py::arg("turbine"), py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z") )
        .def("set_turbine_blade_forces", [](ActuatorFarm& al, uint turbine, arr forcesX, arr forcesY, arr forcesZ)
        { 
            al.setTurbineBladeForces(turbine, static_cast<float *>(forcesX.request().ptr), static_cast<float *>(forcesY.request().ptr), static_cast<float *>(forcesZ.request().ptr)); 
        }, py::arg("turbine"), py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z") )
        .def("calc_blade_forces", &ActuatorFarm::calcBladeForces);
    }
}