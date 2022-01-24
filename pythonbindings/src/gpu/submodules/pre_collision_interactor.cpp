#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/PreCollisionInteractor/PreCollisionInteractor.h>

namespace pre_collision_interactor
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<PreCollisionInteractor, std::shared_ptr<PreCollisionInteractor>>(parentModule, "PreCollisionInteractor");
    }
}