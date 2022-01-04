#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Communication/Communicator.h>

namespace communicator
{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module commModule = parentModule.def_submodule("communicator");

        py::class_<vf::gpu::Communicator>(commModule, "Communicator")
        .def("get_instanz", py::overload_cast<>(&vf::gpu::Communicator::getInstanz))
        .def("get_instanz", py::overload_cast<const int>(&vf::gpu::Communicator::getInstanz))
        .def("get_number_of_process", &vf::gpu::Communicator::getNummberOfProcess)
        .def("get_pid", &vf::gpu::Communicator::getPID);

        return commModule;
    }
}