#ifndef BOOSTSERIALIZATIONCLASSEXPORTHELPER_H
#define BOOSTSERIALIZATIONCLASSEXPORTHELPER_H

#include <D3Q27EsoTwist3DSplittedVector.h>
#include <BCProcessor.h>
#include <BCProcessor.h>
#include <DataSet3D.h>
#include <CompressibleCumulantLBMKernel.h>
#include <IncompressibleCumulantLBMKernel.h>
#include <IncompressibleCumulantWithSpongeLayerLBMKernel.h>
#include <Interactor3D.h>
#include <D3Q27Interactor.h>
#include <Communicator.h>
#include <MPICommunicator.h>
#include <CoProcessor.h>
#include <WriteMacroscopicQuantitiesCoProcessor.h>
#include <ShearStressCoProcessor.h>
#include <AverageValuesCoProcessor.h>
#include <basics/container/CbArray4D.h>
#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <BCAlgorithm.h>
#include <VelocityBCAlgorithm.h>
#include <NonEqDensityBCAlgorithm.h>
#include <EqDensityBCAlgorithm.h>
#include <NoSlipBCAlgorithm.h>
#include <HighViscosityNoSlipBCAlgorithm.h>
#include <SlipBCAlgorithm.h>
#include <NonReflectingDensityBCAlgorithm.h>
#include <NonReflectingVelocityBCAlgorithm.h>

#include <BCAdapter.h>
#include <DensityBCAdapter.h>
#include <NoSlipBCAdapter.h>
#include <SlipBCAdapter.h>
#include <VelocityBCAdapter.h>
#include "ThinWallNoSlipBCAlgorithm.h"
#include "ThinWallBCProcessor.h"
#include "VoidLBMKernel.h"
#include "VoidData3D.h"
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(CompressibleCumulantLBMKernel)
BOOST_CLASS_EXPORT(IncompressibleCumulantLBMKernel)
BOOST_CLASS_EXPORT(IncompressibleCumulantWithSpongeLayerLBMKernel)
BOOST_CLASS_EXPORT(VoidLBMKernel)
BOOST_CLASS_EXPORT(VoidData3D)
BOOST_CLASS_EXPORT(D3Q27EsoTwist3DSplittedVector)
BOOST_CLASS_EXPORT(BCProcessor)
BOOST_CLASS_EXPORT(ThinWallBCProcessor)
BOOST_CLASS_EXPORT(DataSet3D)
BOOST_CLASS_EXPORT(Interactor3D)
BOOST_CLASS_EXPORT(D3Q27Interactor)
BOOST_CLASS_EXPORT(CoProcessor)
BOOST_CLASS_EXPORT(WriteMacroscopicQuantitiesCoProcessor)
BOOST_CLASS_EXPORT(ShearStressCoProcessor)
BOOST_CLASS_EXPORT(AverageValuesCoProcessor)
BOOST_CLASS_EXPORT(WbWriterVtkXmlASCII)
BOOST_CLASS_EXPORT(WbWriterVtkXmlBinary)

//BOOST_CLASS_EXPORT(BCAlgorithm)
//BOOST_CLASS_EXPORT(VelocityBCAlgorithm)
//BOOST_CLASS_EXPORT(NonEqDensityBCAlgorithm)
//BOOST_CLASS_EXPORT(EqDensityBCAlgorithm)
//BOOST_CLASS_EXPORT(NoSlipBCAlgorithm)
//BOOST_CLASS_EXPORT(HighViscosityNoSlipBCAlgorithm)
//BOOST_CLASS_EXPORT(SlipBCAlgorithm)
//BOOST_CLASS_EXPORT(NonReflectingDensityBCAlgorithm)
//BOOST_CLASS_EXPORT(NonReflectingVelocityBCAlgorithm)
//BOOST_CLASS_EXPORT(ThinWallNoSlipBCAlgorithm)

BOOST_CLASS_EXPORT(BCAdapter)
BOOST_CLASS_EXPORT(DensityBCAdapter)
BOOST_CLASS_EXPORT(NoSlipBCAdapter)
BOOST_CLASS_EXPORT(SlipBCAdapter)
BOOST_CLASS_EXPORT(VelocityBCAdapter)
#endif
