#include "TurbulenceModelManager.h"
#include "Parameter/Parameter.h"
#include <memory>


void TurbulenceModelManager::runTurbulenceModelKernel(int level) const
{
    if (this->turbulenceModelKernel.has_value())
        this->turbulenceModelKernel.value()(para.get(), level);
}

void TurbulenceModelManager::runTurbulenceModelADKernel(int level) const
{
    if (this->turbulenceModelADKernel.has_value())
        this->turbulenceModelADKernel.value()(para.get(), level);
}
