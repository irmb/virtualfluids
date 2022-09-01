#include "TurbulenceModelFactory.h"
#include "GPU/TurbulentViscosity.h"
#include "Parameter/Parameter.h"
#include <variant>

void TurbulenceModelFactory::setTurbulenceModel(const TurbulenceModel _turbulenceModel)
{
    this->turbulenceModel = _turbulenceModel;
    para->setTurbulenceModel(_turbulenceModel);
    if(this->turbulenceModel != TurbulenceModel::None) para->setUseTurbulentViscosity(true);

    switch (this->turbulenceModel) {
        case TurbulenceModel::AMD:
            this->turbulenceModelKernel = calcTurbulentViscosityAMD;
            break;
        default:
            this->turbulenceModelKernel = nullptr;
    }
}

void TurbulenceModelFactory::setModelConstant(const real modelConstant)
{
    para->setSGSConstant(modelConstant);
}

void TurbulenceModelFactory::runTurbulenceModelKernel(const int level) const
{
    if(this->turbulenceModelKernel) this->turbulenceModelKernel(para.get(), level);
}
