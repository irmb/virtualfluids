#include "GridScalingFactory.h"
#include "Scaling.cuh"

void GridScalingFactory::setScalingFactory(const GridScalingFactory::GridScaling gridScalingType)
{
    this->gridScalingType = gridScalingType;
}

gridScaling GridScalingFactory::getGridScalingFC(bool hasTurbulentViscosity) const
{
    // for descriptions of the scaling types refer to the header
    switch (gridScalingType) {
        case GridScaling::ScaleRhoSq:
            return ScaleFC_RhoSq_comp_27;
            break;
        case GridScaling::ScaleCompressible:
            if(hasTurbulentViscosity)   return ScaleFC_compressible<true>;
            else                        return ScaleFC_compressible<false>;
            break;
        default:
            return nullptr;
    }
}

gridScaling GridScalingFactory::getGridScalingCF(bool hasTurbulentViscosity) const
{
    // for descriptions of the scaling types refer to the header
    switch (gridScalingType) {
        case GridScaling::ScaleRhoSq:
            return ScaleCF_RhoSq_comp_27;
            break;
        case GridScaling::ScaleCompressible:
            {
                if(hasTurbulentViscosity)   return ScaleCF_compressible<true>;
                else                        return ScaleCF_compressible<false>;
                break;
            }
        default:
            return nullptr;
    }
}