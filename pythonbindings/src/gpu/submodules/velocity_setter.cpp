#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/VelocitySetter.h>

namespace velocity_setter
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<VelocityReader, std::shared_ptr<VelocityReader>>(parentModule, "VelocityReader")
        .def(py::init < std::string,
                        std::string,
                        real,
                        real,
                        real,
                        real>(),
                        "precursor_file_prefix", 
                        "precursor_file_suffix", 
                        "y_start", 
                        "y_end", 
                        "z_start", 
                        "z_end")
        .def("read_slice", &VelocityReader::readSlice);

        py::class_<VelocitySetter, PrecollisionInteractor, std::shared_ptr<VelocitySetter>>(parentModule, "VelocitySetter")
    }
}