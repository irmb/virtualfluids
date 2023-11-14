#include "PreProcessorFactoryImp.h"

#include "PreProcessor/PreProcessorImp.h"

#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionCompressibleD3Q7/InitAdvectionDiffusionCompressibleD3Q7.h"
#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionCompressible/InitAdvectionDiffusionCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitNavierStokesCompressible/InitNavierStokesCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitK18K20NavierStokesCompressible/InitK18K20NavierStokesCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionIncompressible/InitAdvectionDiffusionIncompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionIncompressibleD3Q7/InitAdvectionDiffusionIncompressibleD3Q7.h"
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
    case InitAdvectionDiffusionIncompressibleD3Q7:
        return InitAdvectionDiffusionIncompressibleD3Q7::getNewInstance(para);
        break;
    case InitAdvectionDiffusionIncompressible:
        return InitAdvectionDiffusionIncompressible::getNewInstance(para);
        break;
    case InitAdvectionDiffusionCompressibleD3Q7:
        return InitAdvectionDiffusionCompressibleD3Q7::getNewInstance(para);
        break;
    case InitAdvectionDiffusionCompressible:
        return InitAdvectionDiffusionCompressible::getNewInstance(para);
        break;
    default:
        break;
    }
    throw  std::runtime_error("PreProcessorFactory does not know the PreProcessorType.");
}
