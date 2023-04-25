#ifndef KERNEL_TYPES_H
#define KERNEL_TYPES_H

namespace vf::CollisionKernel::Compressible {
    static const std::string BGK = "BGKCompSP27";
    static const std::string BGKUnified = "BGKUnified";
    static const std::string BGKPlus = "BGKPlusCompSP27";
    static const std::string MRT = "MRTCompSP27";
    static const std::string Cascade = "CascadeCompSP27";

    static const std::string CumulantClassic = "CumulantCompSP27";

    static const std::string CumulantK15Unified = "CumulantK15Unified";
    static const std::string CumulantK17Unified = "CumulantK17Unified";

    static const std::string CumulantK17Bulk = "CumulantK17BulkComp";
    static const std::string CumulantK17Chim = "CumulantK17CompChim";
    static const std::string CumulantK17 = "CumulantK17";

    static const std::string CumulantAll4SP27 = "CumulantAll4CompSP27";
    static const std::string CumulantK18 = "CumulantK18Comp";
    static const std::string CumulantK20 = "CumulantK20Comp";

    static const std::string CumulantK15 = "CumulantK15Comp";
    static const std::string CumulantK15Bulk = "CumulantK15BulkComp";
    static const std::string CumulantK15Sponge = "CumulantK15SpongeComp";
}

namespace vf::CollisionKernel::Incompressible {
    static const std::string BGK = "BGKIncompSP27";
    static const std::string BGKPlus = "BGKPlusIncompSP27";
    static const std::string MRT = "MRTIncompSP27";
    static const std::string Cascade = "CascadeIncompSP27";

    static const std::string Cumulant1h = "Cumulant1hIncompSP27";
    static const std::string CumulantIsometric = "CumulantIsoIncompSP27";
    static const std::string CumulantK15 = "CumulantK15Incomp";
}

namespace vf::CollisionKernel::PorousMedia {
    static const std::string CumulantOne = "CumulantOneCompSP27";
}

namespace vf::CollisionKernel::Wale {
    static const std::string CumulantK17 = "WaleCumulantK17Comp";
    static const std::string CumulantK17Debug = "WaleCumulantK17DebugComp";
    static const std::string CumulantK15 = "WaleCumulantK15Comp";
    static const std::string CumulantK15SoniMalav = "WaleBySoniMalavCumulantK15Comp";
}

#endif
