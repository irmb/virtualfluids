#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/Probe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace probes
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module probeModule = parentModule.def_submodule("probes");

        py::enum_<PostProcessingVariable>(probeModule, "PostProcessingVariables")
        .value("Means", PostProcessingVariable::Means)
        .value("Variances", PostProcessingVariable::Variances);

        py::class_<Probe, PreCollisionInteractor, std::shared_ptr<Probe>>(probeModule, "Probe")
        .def("add_post_processing_variable", &Probe::addPostProcessingVariable);

        py::class_<PointProbe, Probe, std::shared_ptr<PointProbe>>(probeModule, "PointProbe")
        .def(py::init<
                        const std::string,
                        uint,
                        uint, 
                        uint>(), 
                        "probe_name",
                        "t_start_avg",
                        "t_start_out",
                        "t_out")
        .def("add_probe_points_from_list", &PointProbe::addProbePointsFromList)
        .def("add_probe_points_from_x_normal_plane", &PointProbe::addProbePointsFromXNormalPlane);

        py::class_<PlaneProbe, Probe, std::shared_ptr<PlaneProbe>>(probeModule, "PlaneProbe")
        .def(py::init<
                        const std::string,
                        uint,
                        uint, 
                        uint>(), 
                        "probe_name",
                        "t_start_avg",
                        "t_start_out",
                        "t_out")
        .def("set_probe_plane", &PlaneProbe::setProbePlane);

        return probeModule;
    }
}