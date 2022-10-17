#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <gpu/VirtualFluids_GPU/Parameter/Parameter.h>
#include <basics/config/ConfigurationFile.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace parameter
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<Parameter, std::shared_ptr<Parameter>>(parentModule, "Parameter")
        .def(py::init<
                int,
                int,
                std::optional<const vf::basics::ConfigurationFile*>>(),
                "number_of_processes",
                "my_ID",
                "config_data")
        .def(py::init<int, int>(),
                "number_of_processes",
                "my_ID")
        .def(py::init<const vf::basics::ConfigurationFile*>(), "config_data")
        .def("set_forcing", &Parameter::setForcing)
        .def("set_quadric_limiters", &Parameter::setQuadricLimiters)
        .def("set_diff_on", &Parameter::setDiffOn)
        .def("set_comp_on", &Parameter::setCompOn)
        .def("set_max_level", &Parameter::setMaxLevel)
        .def("set_timestep_end", &Parameter::setTimestepEnd)
        .def("set_timestep_out", &Parameter::setTimestepOut)
        .def("set_timestep_start_out", &Parameter::setTimestepStartOut)
        .def("set_timestep_of_coarse_level", &Parameter::setTimestepOfCoarseLevel)
        .def("set_calc_turbulence_intensity", &Parameter::setCalcTurbulenceIntensity)
        .def("set_output_path", &Parameter::setOutputPath)
        .def("set_output_prefix", &Parameter::setOutputPrefix)
        .def("set_print_files", &Parameter::setOutflowPressureCorrectionFactor)
        .def("set_print_files", &Parameter::setPrintFiles)
        .def("set_temperature_init", &Parameter::setTemperatureInit)
        .def("set_temperature_BC", &Parameter::setTemperatureBC)
        .def("set_viscosity_LB", &Parameter::setViscosityLB)
        .def("set_velocity_LB", &Parameter::setVelocityLB)
        .def("set_viscosity_ratio", &Parameter::setViscosityRatio)
        .def("set_velocity_ratio", &Parameter::setVelocityRatio)
        .def("set_density_ratio", &Parameter::setDensityRatio)
        .def("set_devices", &Parameter::setDevices)
        .def("set_is_body_force", &Parameter::setIsBodyForce)
        .def("set_use_streams", &Parameter::setUseStreams)
        .def("set_main_kernel", &Parameter::setMainKernel)
        .def("set_AD_kernel", &Parameter::setADKernel)
        .def("set_has_wall_monitor", &Parameter::setHasWallModelMonitor)
        .def("set_outflow_pressure_correction_factor", &Parameter::setOutflowPressureCorrectionFactor)
        .def("set_initial_condition", [](Parameter &para, std::function<std::vector<float>(real, real, real)> &init_func)
        {
            para.setInitialCondition([init_func](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
            {   
                std::vector<float> values = init_func(coordX, coordY, coordZ);
                rho = values[0];
                vx = values[1];
                vy = values[2];
                vz = values[3];
            });
        })
        .def("set_initial_condition_uniform", [](Parameter &para, real velocity_x, real velocity_y, real velocity_z)
        {
            para.setInitialCondition([velocity_x, velocity_y, velocity_z](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) // must capture values explicitly!
            {
                rho = c0o1;
                vx = velocity_x;
                vy = velocity_y;
                vz = velocity_z;
            });
        })
        .def("set_initial_condition_log_law", [](Parameter &para, real u_star, real z0, real velocityRatio)
        {
            para.setInitialCondition(
                [u_star, z0, velocityRatio](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
                {
                    coordZ = coordZ > c0o1 ? coordZ : c0o1;

                    rho = c0o1;
                    vx  = u_star/c4o10 * log(coordZ/z0+c1o1) / velocityRatio;
                    vy = c0o1;
                    vz = c0o1;
                }
            );
        })
        .def("set_initial_condition_perturbed_log_law", [](Parameter &para, real u_star, real z0, real L_x, real L_z, real H, real velocityRatio)
        {
            para.setInitialCondition(
                [u_star, z0, L_x, L_z, H, velocityRatio](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
                {
                    coordZ = coordZ > c0o1 ? coordZ : c0o1;
                    rho = c0o1;
                    vx  = (u_star/c4o10 * log(coordZ/z0+c1o1) + c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1)) / velocityRatio; 
                    vy  = c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1) / velocityRatio; 
                    vz  = c8o1*u_star/c4o10*(sin(cPi*c8o1*coordY/H)*sin(cPi*c8o1*coordZ/H)+sin(cPi*c8o1*coordX/L_x))/(pow(c1o2*L_z-coordZ, c2o1)+c1o1) / velocityRatio;
                }
            );
        })
        .def("set_has_wall_model_monitor", &Parameter::setHasWallModelMonitor)
        .def("add_actuator", &Parameter::addActuator)
        .def("add_probe", &Parameter::addProbe)
        .def("get_output_path", &Parameter::getOutputPath)
        .def("get_output_prefix", &Parameter::getOutputPrefix)
        .def("get_velocity", &Parameter::getVelocity)
        .def("get_viscosity", &Parameter::getViscosity)
        .def("get_velocity_ratio", &Parameter::getVelocityRatio)
        .def("get_viscosity_ratio", &Parameter::getViscosityRatio)
        .def("get_density_ratio", &Parameter::getDensityRatio)
        .def("get_force_ratio", &Parameter::getForceRatio)
        .def("get_SGS_constant", &Parameter::getSGSConstant)
        .def("get_is_body_force", &Parameter::getIsBodyForce)
        ;

    }
}