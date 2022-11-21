#include <pybind11/pybind11.h>
#include "basics/basics.cpp"
#include "lbm/lbm.cpp"
#include "logger/logger.cpp"

#ifdef VF_GPU_PYTHONBINDINGS
#include "gpu/gpu.cpp"
#endif
#ifdef VF_CPU_PYTHONBINDINGS
#include "cpu/cpu.cpp"
#endif


namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(bindings, m)
    {
        py::add_ostream_redirect(m, "ostream_redirect");
        basics::makeModule(m);
        lbm::makeModule(m);
        logging::makeModule(m);
#ifdef VF_GPU_PYTHONBINDINGS
        gpu::makeModule(m);
#endif
#ifdef VF_CPU_PYTHONBINDINGS
        cpu::makeModule(m);
#endif
    }
}