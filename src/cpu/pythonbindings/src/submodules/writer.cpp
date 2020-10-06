#include <pybind11/pybind11.h>
#include <simulationconfig/WriterConfig.h>

namespace py = pybind11;

void makeWriterModule(py::module &parentModule)
{
    py::module writerModule = parentModule.def_submodule("writer");

    py::enum_<WriterType>(writerModule, "WriterType")
            .value("ASCII", WriterType::ASCII)
            .value("BINARY", WriterType::BINARY);

    py::class_<WriterConfig>(writerModule, "Writer")
            .def(py::init())
            .def_readwrite("output_path", &WriterConfig::outputPath)
            .def_readwrite("type", &WriterConfig::writerType);
}