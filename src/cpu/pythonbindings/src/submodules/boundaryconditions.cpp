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
using namespace py::literals;

template<class Adapter>
using py_bc_class = py::class_<Adapter, BCAdapter, std::shared_ptr<Adapter>>;

template<class Adapter, class Algorithm, typename ...Args>
std::shared_ptr<Adapter> create(Args... args)
{
    auto adapter = std::make_shared<Adapter>(args...);
    adapter->setBcAlgorithm(std::make_shared<Algorithm>());
    return adapter;
}

template<class Algorithm>
void add_constructors_to_velocity_bc(py_bc_class<VelocityBCAdapter> &cls)
{
    auto firstConstructor = &create<VelocityBCAdapter, Algorithm, bool &, bool &, bool &, mu::Parser &, double &, double &>;
    auto secondConstructor = &create<VelocityBCAdapter, Algorithm,
            bool &, bool &, bool &, mu::Parser &, mu::Parser &, mu::Parser &, double &, double &>;
    auto thirdConstructor = &create<VelocityBCAdapter, Algorithm,
            double &, double &, double &, double &, double &, double &, double &, double &, double &>;

    cls.def(py::init(&create<VelocityBCAdapter, Algorithm>))
            .def(py::init(firstConstructor),
                 "vx1"_a, "vx2"_a, "vx3"_a, "function"_a, "start_time"_a, "end_time"_a)
            .def(py::init(secondConstructor),
                 "vx1"_a, "vx2"_a, "vx3"_a,
                 "function_vx1"_a, "function_vx2"_a, "function_vx2"_a,
                 "start_time"_a, "end_time"_a)
            .def(py::init(thirdConstructor),
                 "vx1"_a, "vx1_start_time"_a, "vx1_end_time"_a,
                 "vx2"_a, "vx2_start_time"_a, "vx2_end_time"_a,
                 "vx3"_a, "vx3_start_time"_a, "vx3_end_time"_a);
}

void makeBoundaryConditionsModule(py::module_ &parentModule)
{
    py::module_ bcModule = parentModule.def_submodule("boundaryconditions");

    py::class_<BCAdapter, std::shared_ptr<BCAdapter>>(bcModule, "BCAdapter");

    py_bc_class<NoSlipBCAdapter>(bcModule, "NoSlipBoundaryCondition")
            .def(py::init(&create<NoSlipBCAdapter, NoSlipBCAlgorithm>));

    auto velocityBoundaryCondition = py_bc_class<VelocityBCAdapter>(bcModule, "VelocityBoundaryCondition");
    add_constructors_to_velocity_bc<VelocityBCAlgorithm>(velocityBoundaryCondition);

    py_bc_class<DensityBCAdapter>(bcModule, "DensityBoundaryCondition")
            .def(py::init(&create<DensityBCAdapter, NonReflectingOutflowBCAlgorithm>));
}

