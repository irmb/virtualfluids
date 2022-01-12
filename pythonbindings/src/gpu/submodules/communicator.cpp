#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Communication/Communicator.h>

namespace communicator
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<vf::gpu::Communicator>(parentModule, "Communicator")
        .def("get_instanz", py::overload_cast<>(&vf::gpu::Communicator::getInstanz), py::return_value_policy::reference)
        // .def("get_instanz", py::overload_cast<const int>(&vf::gpu::Communicator::getInstanz), py::return_value_policy::reference)
        .def("get_number_of_process", &vf::gpu::Communicator::getNummberOfProcess)
        .def("get_pid", &vf::gpu::Communicator::getPID);
    }
}