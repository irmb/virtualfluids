#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Communication/Communicator.h>

namespace communicator
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<vf::gpu::Communicator, std::unique_ptr<vf::gpu::Communicator, py::nodelete>>(parentModule, "Communicator")
        .def_static("get_instance", &vf::gpu::Communicator::getInstance, py::return_value_policy::reference)
        .def("get_number_of_process", &vf::gpu::Communicator::getNummberOfProcess)
        .def("get_pid", &vf::gpu::Communicator::getPID);
    }
}