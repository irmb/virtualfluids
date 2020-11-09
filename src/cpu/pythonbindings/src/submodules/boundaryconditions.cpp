#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <BoundaryConditions/DensityBCAdapter.h>
#include <BoundaryConditions/NonReflectingOutflowBCAlgorithm.h>
#include <BoundaryConditions/BCAdapter.h>
#include <BoundaryConditions/NoSlipBCAdapter.h>
#include <BoundaryConditions/VelocityBCAdapter.h>
#include <BoundaryConditions/NoSlipBCAlgorithm.h>
#include <BoundaryConditions/VelocityBCAlgorithm.h>

namespace py = pybind11;

void makeBoundaryConditionsModule(py::module_ &parentModule)
{
    py::module bcModule = parentModule.def_submodule("boundaryconditions");

    py::class_<BCAdapter, std::shared_ptr<BCAdapter>>(bcModule, "BCAdapter")
            .def_property("algorithm", &BCAdapter::getAlgorithm, &BCAdapter::setBcAlgorithm);

    py::class_<NoSlipBCAdapter, BCAdapter, std::shared_ptr<NoSlipBCAdapter>>(bcModule, "NoSlipBCAdapter")
            .def(py::init());

    py::class_<VelocityBCAdapter, BCAdapter, std::shared_ptr<VelocityBCAdapter>>(bcModule, "VelocityBCAdapter")
            .def(py::init())
            .def(py::init<bool &, bool &, bool &, mu::Parser &, double &, double &>())
            .def(py::init<bool &, bool &, bool &, mu::Parser &, mu::Parser &, mu::Parser &, double &, double &>())
            .def(py::init<bool &, bool &, bool &, std::string &, double &, double &>())
            .def(py::init<BCFunction &, bool, bool, bool>())
            .def(py::init<BCFunction &, BCFunction &, BCFunction &>())
            .def(py::init<std::vector<BCFunction> &, std::vector<BCFunction> &, std::vector<BCFunction> &>())
            .def(py::init<double &, double &, double &, double &, double &, double &, double &, double &, double &>())
            .def(py::init<std::string &, double &, double &, std::string &, double &, double &, std::string &, double &, double &>());

    py::class_<DensityBCAdapter, BCAdapter, std::shared_ptr<DensityBCAdapter>>(bcModule, "DensityBCAdapter")
            .def(py::init());


    py::class_<BCAlgorithm, std::shared_ptr<BCAlgorithm>>(bcModule, "BCAlgorithm");

    py::class_<NoSlipBCAlgorithm, BCAlgorithm, std::shared_ptr<NoSlipBCAlgorithm>>(bcModule, "NoSlipBCAlgorithm")
            .def(py::init());

    py::class_<VelocityBCAlgorithm, BCAlgorithm, std::shared_ptr<VelocityBCAlgorithm>>(bcModule, "VelocityBCAlgorithm")
            .def(py::init());

    py::class_<NonReflectingOutflowBCAlgorithm, BCAlgorithm, std::shared_ptr<NonReflectingOutflowBCAlgorithm>>(bcModule,
                                                                                                               "NonReflectingOutflowBCAlgorithm")
            .def(py::init());
}

