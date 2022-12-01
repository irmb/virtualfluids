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
        .def("add_statistic", &Probe::addStatistic, py::arg("variable"))
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
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"))
        .def("add_probe_points_from_list", &PointProbe::addProbePointsFromList, py::arg("point_coords_x"), py::arg("point_coords_y"), py::arg("point_coords_z"))
        .def("add_probe_points_from_x_normal_plane", &PointProbe::addProbePointsFromXNormalPlane, py::arg("pos_x"), py::arg("pos0_y"), py::arg("pos0_z"), py::arg("pos1_y"), py::arg("pos1_z"), py::arg("n_y"), py::arg("n_z"));

        py::class_<PlaneProbe, Probe, std::shared_ptr<PlaneProbe>>(probeModule, "PlaneProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint>(), 
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"))
        .def("set_probe_plane", &PlaneProbe::setProbePlane, py::arg("pos_x"), py::arg("pos_y"), py::arg("pos_z"), py::arg("delta_x"), py::arg("delta_y"), py::arg("delta_z"));

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
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("plane_normal"));


        py::class_<WallModelProbe, Probe, std::shared_ptr<WallModelProbe>>(probeModule, "WallModelProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        uint>(), 
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"))
        .def("set_force_output_to_stress", &WallModelProbe::setForceOutputToStress, py::arg("output_stress"))
        .def("set_evaluate_pressure_gradient", &WallModelProbe::setEvaluatePressureGradient, py::arg("eval_press_grad"));

        return probeModule;
    }
}