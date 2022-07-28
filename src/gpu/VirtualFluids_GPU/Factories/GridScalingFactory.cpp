#include "GridScalingFactory.h"
#include "GPU/GPU_Interface.h"

void GridScalingFactory::setScalingFactory(const GridScalingFactory::GridScaling gridScalingType)
{
    this->gridScaling = gridScalingType;
}

gridScalingFC GridScalingFactory::getGridScalingFC() const
{
    // for descriptions of the scaling types refer to the header
    switch (gridScaling) {
        case GridScaling::ScaleRhoSq:
            return ScaleFC_RhoSq_comp_27;
            break;
        case GridScaling::ScaleK17:
            return ScaleFC_K17_redesigned;
            break;
        default:
            return nullptr;
    }
}

// gridScalingCF GridScalingFactory::getGridScalingCF() const
// {
//     // for descriptions of the scaling types refer to the header
//     switch (gridScaling) {
//         case GridScaling::ScaleRhoSq:
//             return ScaleCF_RhoSq_comp_27;
//             break;
//         case GridScaling::ScaleK17:
//             return ScaleCF_K17_redesigned;
//             break;
//         default:
//             return nullptr;
//     }
// }