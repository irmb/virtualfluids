#include <pybind11/pybind11.h>
#include <gpu/GridGenerator/VelocitySetter/VelocitySetter.h>

namespace velocity_setter
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<FileType>(parentModule, "FileType")
        .value("VTK", FileType::VTK);

        parentModule.def("create_file_collection", &createFileCollection);

        py::class_<VelocityFileCollection, std::shared_ptr<VelocityFileCollection>>(parentModule, "VelocityFileCollection");

        py::class_<VTKFileCollection, VelocityFileCollection, std::shared_ptr<VTKFileCollection>>(parentModule, "VTKFileCollection")
        .def(py::init <std::string>(), "prefix");
    }
}