#include "PreProcessorFactoryImp.h"

#include "PreProcessor/PreProcessorImp.h"

#include "PreProcessor/PreProcessorStrategy/InitCompAD7/InitCompAD7.h"
#include "PreProcessor/PreProcessorStrategy/InitCompAD27/InitCompAD27.h"
#include "PreProcessor/PreProcessorStrategy/InitNavierStokesCompressible/InitNavierStokesCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitK18K20NavierStokesCompressible/InitK18K20NavierStokesCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitIncompAD27/InitIncompAD27.h"
#include "PreProcessor/PreProcessorStrategy/InitIncompAD7/InitIncompAD7.h"
#include "PreProcessor/PreProcessorStrategy/InitNavierStokesIncompressible/InitNavierStokesIncompressible.h"


std::shared_ptr<PreProcessor> PreProcessorFactoryImp::makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para)
{
    std::shared_ptr<PreProcessorImp> prePro = PreProcessorImp::getNewInstance();

    for (std::size_t i = 0; i < preProcessorTypes.size(); i++)
        prePro->addStrategy(makePreProcessorStrategy(preProcessorTypes.at(i), para));

    return prePro;
}

std::shared_ptr<PreProcessorStrategy> PreProcessorFactoryImp::makePreProcessorStrategy(PreProcessorType preProcessorType, std::shared_ptr<Parameter> para)
{
    switch (preProcessorType)
    {
    case InitNavierStokesIncompressible:
        return InitNavierStokesIncompressible::getNewInstance(para);
        break;
    case InitNavierStokesCompressible:
        return InitNavierStokesCompressible::getNewInstance(para);
        break;
    case InitK18K20NavierStokesCompressible:
        return InitK18K20NavierStokesCompressible::getNewInstance(para);
        break;
    case InitIncompAD7:
        return InitIncompAD7::getNewInstance(para);
        break;
    case InitIncompAD27:
        return InitIncompAD27::getNewInstance(para);
        break;
    case InitCompAD7:
        return InitCompAD7::getNewInstance(para);
        break;
    case InitCompAD27:
        return InitCompAD27::getNewInstance(para);
        break;
    default:
        break;
    }
    throw  std::runtime_error("PreProcessorFactory does not know the PreProcessorType.");
}
