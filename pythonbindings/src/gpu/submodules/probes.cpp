#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/Probe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/WallModelProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/Probes/PlanarAverageProbe.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace probes
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module probeModule = parentModule.def_submodule("probes");

        py::enum_<Statistic>(probeModule, "Statistic")
        .value("Instantaneous", Statistic::Instantaneous)
        .value("Means", Statistic::Means)
        .value("Variances", Statistic::Variances)
        .value("SpatialMeans", Statistic::SpatialMeans)
        .value("SpatioTemporalMeans", Statistic::SpatioTemporalMeans)
        .value("SpatialCovariances", Statistic::SpatialCovariances)
        .value("SpatioTemporalCovariances", Statistic::SpatioTemporalCovariances)
        .value("SpatialSkewness", Statistic::SpatialSkewness)
        .value("SpatioTemporalSkewness", Statistic::SpatioTemporalSkewness)
        .value("SpatialFlatness", Statistic::SpatialFlatness)
        .value("SpatioTemporalFlatness", Statistic::SpatioTemporalFlatness);

        py::class_<Probe, PreCollisionInteractor, std::shared_ptr<Probe>>(probeModule, "Probe")
        .def("add_statistic", &Probe::addStatistic)
        .def("set_file_name_to_n_out", &Probe::setFileNameToNOut)
        .def("add_all_available_statistics", &Probe::addAllAvailableStatistics);

        py::class_<PointProbe, Probe, std::shared_ptr<PointProbe>>(probeModule, "PointProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint>(), 
                        "probe_name",
                        "output_path"
                        "t_start_avg",
                        "t_avg",
                        "t_start_out",
                        "t_out")
        .def("add_probe_points_from_list", &PointProbe::addProbePointsFromList)
        .def("add_probe_points_from_x_normal_plane", &PointProbe::addProbePointsFromXNormalPlane);

        py::class_<PlaneProbe, Probe, std::shared_ptr<PlaneProbe>>(probeModule, "PlaneProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint>(), 
                        "probe_name",
                        "output_path"
                        "t_start_avg",
                        "t_avg",
                        "t_start_out",
                        "t_out")
        .def("set_probe_plane", &PlaneProbe::setProbePlane);

        py::class_<PlanarAverageProbe, Probe, std::shared_ptr<PlanarAverageProbe>>(probeModule, "PlanarAverageProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint,
                        uint,
                        uint,
                        uint,
                        char>(),
                        "probe_name",
                        "output_path",
                        "t_start_avg",
                        "t_start_tmp_avg",
                        "t_avg",
                        "t_start_out",
                        "t_out",
                        "plane_normal");


        py::class_<WallModelProbe, Probe, std::shared_ptr<WallModelProbe>>(probeModule, "WallModelProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        uint>(), 
                        "probe_name",
                        "output_path"
                        "t_start_avg",
                        "t_start_tmp_avg",
                        "t_avg",
                        "t_start_out",
                        "t_out")
        .def("set_force_output_to_stress", &WallModelProbe::setForceOutputToStress)
        .def("set_evaluate_pressure_gradient", &WallModelProbe::setEvaluatePressureGradient);

        return probeModule;
    }
}