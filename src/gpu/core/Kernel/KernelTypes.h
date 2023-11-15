#ifndef KERNEL_TYPES_H
#define KERNEL_TYPES_H

namespace vf::collisionKernel::compressible {
    static const std::string BGK = "BGKCompSP27";
    static const std::string BGKPlus = "BGKPlusCompSP27";
    static const std::string K17CompressibleNavierStokes = "K17CompressibleNavierStokes";
    static const std::string K15CompressibleNavierStokes = "K15CompressibleNavierStokes";
    }

namespace vf::collisionKernel::incompressible {
    static const std::string BGK = "BGKIncompSP27";
    static const std::string BGKPlus = "BGKPlusIncompSP27";
    static const std::string CumulantK15 = "CumulantK15Incomp";
}

#endif
