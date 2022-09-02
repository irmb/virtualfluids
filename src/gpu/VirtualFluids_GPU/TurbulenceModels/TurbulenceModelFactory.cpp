#include "TurbulenceModelFactory.h"
#include "GPU/TurbulentViscosity.h"
#include "Parameter/Parameter.h"
#include <logger/Logger.h>

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

void TurbulenceModelFactory::readConfigFile(const vf::basics::ConfigurationFile &configData)
{
    if (configData.contains("TurbulenceModel"))
    {
        std::string config = configData.getValue<std::string>("TurbulenceModel");
        
        if      (config == "Smagorinsky") this->setTurbulenceModel( TurbulenceModel::Smagorinsky );
        else if (config == "AMD")         this->setTurbulenceModel( TurbulenceModel::AMD );
        else if (config == "QR" )         this->setTurbulenceModel( TurbulenceModel::QR );
        else if (config == "None")        this->setTurbulenceModel( TurbulenceModel::None );
        else    std::runtime_error("TurbulenceModelFactory: Invalid turbulence model!");

        VF_LOG_INFO("Turbulence model: {}", config);
        
    }

    if (configData.contains("SGSconstant"))
    {
        para->setSGSConstant(configData.getValue<real>("SGSconstant"));

        VF_LOG_INFO("SGS constant: {}", para->getSGSConstant() );
    }
}

void TurbulenceModelFactory::runTurbulenceModelKernel(const int level) const
{
    if(this->turbulenceModelKernel) this->turbulenceModelKernel(para.get(), level);
}
