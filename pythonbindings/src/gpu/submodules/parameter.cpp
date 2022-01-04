#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Parameter/Parameter.h>
#include <basics/config/ConfigurationFile.h>

namespace parameter
{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module paramModule = parentModule.def_submodule("parameter");

        py::class_<Parameter>(paramModule, "Parameter")
        .def(py::init<
                const vf::basics::ConfigurationFile&, 
                int,
                int
                >(),
                "config_data",
                "number_of_processes",
                "my_ID")
        .def("set_forcing", &Parameter::setForcing)
        .def("set_diff_on", &Parameter::setDiffOn)
        .def("set_comp_on", &Parameter::setCompOn)
        .def("set_max_level", &Parameter::setMaxLevel)
        .def("set_t_end", &Parameter::setTEnd)
        .def("set_t_out", &Parameter::setTOut)
        .def("set_t_start_out", &Parameter::setTStartOut)
        .def("set_timestep_of_coarse_level", &Parameter::setTimestepOfCoarseLevel)
        .def("set_output_path", &Parameter::setOutputPath)
        .def("set_output_prefix", &Parameter::setOutputPrefix)
        .def("set_print_files", &Parameter::setPrintFiles)
        .def("set_temperature_init", &Parameter::setTemperatureInit)
        .def("set_temperature_BC", &Parameter::setTemperatureBC)
        .def("set_temperature_init", &Parameter::setTemperatureInit)
        .def("set_viscosity", &Parameter::setViscosity)
        .def("set_velocity", &Parameter::setVelocity)
        .def("set_viscosity_ratio", &Parameter::setViscosityRatio)
        .def("set_velocity_ratio", &Parameter::setVelocityRatio)
        .def("set_devices", &Parameter::setDevices)
        .def("set_is_body_force", &Parameter::setIsBodyForce)
        .def("set_main_kernel", &Parameter::setMainKernel)
        .def("set_AD_kernel", &Parameter::setADKernel)
        .def("add_actuator", &Parameter::addActuator)
        .def("add_probe", &Parameter::addProbe);

        return paramModule;
    }
}