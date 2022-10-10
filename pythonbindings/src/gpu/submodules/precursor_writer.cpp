#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PrecursorWriter.h>

namespace precursor_writer
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<OutputVariable>(parentModule, "OutputVariable")
        .value("Velocities", OutputVariable::Velocities)
        .value("Distributions", OutputVariable::Distributions);

        py::class_<PrecursorWriter, PreCollisionInteractor, std::shared_ptr<PrecursorWriter>>(parentModule, "PrecursorWriter")
        .def(py::init < std::string,
                        std::string,
                        real,
                        real, real,
                        real, real,
                        uint, uint, 
                        OutputVariable, 
                        uint>(),
                        "filename"
                        "output_path", 
                        "x_pos",
                        "y_min", "y_max",
                        "z_min", "z_max",
                        "t_start_out", "t_save", 
                        "output_variable", 
                        "max_timesteps_per_file");
    }
}