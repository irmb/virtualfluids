#ifndef KERNEL_TYPES_H
#define KERNEL_TYPES_H

namespace vf::collisionKernel::compressible {
    static const std::string BGK = "BGKCompSP27";
    static const std::string BGKUnified = "BGKUnified";
    static const std::string BGKPlus = "BGKPlusCompSP27";
    static const std::string MRT = "MRTCompSP27";
    static const std::string Cascade = "CascadeCompSP27";

    static const std::string CumulantClassic = "CumulantCompSP27";

    static const std::string CumulantK15Unified = "CumulantK15Unified";
    static const std::string K17CompressibleNavierStokesUnified = "K17CompressibleNavierStokesUnified";

    static const std::string K17CompressibleNavierStokes = "K17CompressibleNavierStokes";
    static const std::string K17CompressibleNavierStokesBulkViscosity = "K17CompressibleNavierStokesBulkViscosity";
    static const std::string K17CompressibleNavierStokesChimeraLegacy = "K17CompressibleNavierStokesChimeraLegacy";

    static const std::string CumulantAll4SP27 = "CumulantAll4CompSP27";
    static const std::string CumulantK18 = "CumulantK18Comp";
    static const std::string CumulantK20 = "CumulantK20Comp";

    static const std::string K15CompressibleNavierStokes = "K15CompressibleNavierStokes";
    static const std::string K15CompressibleNavierStokesBulk = "K15CompressibleNavierStokesBulk";
    static const std::string K15CompressibleNavierStokesSponge = "K15CompressibleNavierStokesSponge";
    }

namespace vf::collisionKernel::incompressible {
    static const std::string BGK = "BGKIncompSP27";
    static const std::string BGKPlus = "BGKPlusIncompSP27";
    static const std::string MRT = "MRTIncompSP27";
    static const std::string Cascade = "CascadeIncompSP27";

    static const std::string Cumulant1h = "Cumulant1hIncompSP27";
    static const std::string CumulantIsometric = "CumulantIsoIncompSP27";
    static const std::string CumulantK15 = "CumulantK15Incomp";
}

namespace vf::collisionKernel::porousMedia {
    static const std::string CumulantOne = "CumulantOneCompSP27";
}

namespace vf::collisionKernel::wale {
    static const std::string CumulantK17 = "WaleCumulantK17Comp";
    static const std::string CumulantK17Debug = "WaleCumulantK17DebugComp";
    static const std::string CumulantK15 = "WaleCumulantK15Comp";
    static const std::string CumulantK15SoniMalav = "WaleBySoniMalavCumulantK15Comp";
}

#endif