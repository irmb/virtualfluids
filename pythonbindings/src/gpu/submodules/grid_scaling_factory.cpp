#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Factories/GridScalingFactory.h>

namespace grid_scaling_factory
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        
        py::class_<GridScalingFactory, std::shared_ptr<GridScalingFactory>>(parentModule, "GridScalingFactory")
        .def(py::init<>())
        .def("set_scaling_factory", &GridScalingFactory::setScalingFactory, py::arg("scaling_type"));

        py::enum_<GridScalingFactory::GridScaling>(parentModule, "GridScaling")
        .value("ScaleCompressible", GridScalingFactory::GridScaling::ScaleCompressible)
        .value("ScaleRhoSq", GridScalingFactory::GridScaling::ScaleRhoSq)
        .value("NotSpecified", GridScalingFactory::GridScaling::NotSpecified);
    }
}