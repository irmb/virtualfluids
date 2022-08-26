#include "TurbulenceModelFactory.h"
#include "GPU/TurbulentViscosity.h"
#include "Parameter/Parameter.h"
#include <variant>

void TurbulenceModelFactory::setTurbulenceModel(const TurbulenceModelFactory::TurbulenceModel _turbulenceModel)
{
    this->turbulenceModel = _turbulenceModel;
}

TurbulenceModelFactory::TurbulenceModel TurbulenceModelFactory::getTurbulenceModel()
{
    return this->turbulenceModel;
}

TurbulenceModelKernel TurbulenceModelFactory::getTurbulenceModelKernel() const
{
    switch (this->turbulenceModel) {
        case TurbulenceModel::AMD:
            return calcTurbulentViscosityAMD;
            break;
        default:
            return nullptr;
    }
}


