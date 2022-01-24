#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Output/FileWriter.h>
#include <gpu/VirtualFluids_GPU/Output/DataWriter.h>


namespace file_writer
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<DataWriter, std::shared_ptr<DataWriter>>(parentModule, "_DataWriter");

        py::class_<FileWriter, DataWriter, std::shared_ptr<FileWriter>>(parentModule, "FileWriter")
        .def(py::init<>());
    }
}