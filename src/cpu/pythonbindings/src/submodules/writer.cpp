#include <pybind11/pybind11.h>
#include <simulationconfig/WriterConfiguration.h>

namespace py = pybind11;

void makeWriterModule(py::module_ &parentModule)
{
    py::module writerModule = parentModule.def_submodule("writer");

    py::enum_<OutputFormat>(writerModule, "OutputFormat")
            .value("ASCII", OutputFormat::ASCII)
            .value("BINARY", OutputFormat::BINARY);

    py::class_<WriterConfiguration>(writerModule, "Writer")
            .def(py::init())
            .def_readwrite("output_path", &WriterConfiguration::outputPath)
            .def_readwrite("output_format", &WriterConfiguration::outputFormat);
}