#ifndef BOOSTSERIALIZATIONCLASSEXPORTHELPER_H
#define BOOSTSERIALIZATIONCLASSEXPORTHELPER_H

#include <LBMKernelETD3Q27.h>
#include <D3Q27EsoTwist3DSplittedVector.h>
#include <BCProcessor.h>
#include <D3Q27ETBCProcessor.h>
#include <DataSet3D.h>
#include <LBMKernelETD3Q27CCLB.h>
#include <LBMKernelETD3Q27CCLBWithSpongeLayer.h>
#include <Interactor3D.h>
#include <D3Q27Interactor.h>
#include <Communicator.h>
#include <MPICommunicator.h>
#include <CoProcessor.h>
#include <MacroscopicQuantitiesCoProcessor.h>
#include <ShearStressCoProcessor.h>
#include <AverageValuesCoProcessor.h>
#include <basics/container/CbArray4D.h>
#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <BoundaryCondition.h>
#include <VelocityBoundaryCondition.h>
#include <NonEqDensityBoundaryCondition.h>
#include <EqDensityBoundaryCondition.h>
#include <NoSlipBoundaryCondition.h>
#include <HighViscosityNoSlipBoundaryCondition.h>
#include <SlipBoundaryCondition.h>
#include <NonReflectingDensityBoundaryCondition.h>
#include <NonReflectingVelocityBoundaryCondition.h>

#include <D3Q27BoundaryConditionAdapter.h>
#include <D3Q27DensityBCAdapter.h>
#include <D3Q27NoSlipBCAdapter.h>
#include <D3Q27SlipBCAdapter.h>
#include <D3Q27VelocityBCAdapter.h>
#include "ThinWallNoSlipBoundaryCondition.h"
#include "D3Q27ETForThinWallBCProcessor.h"

#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(LBMKernelETD3Q27)
BOOST_CLASS_EXPORT(LBMKernelETD3Q27CCLB)
BOOST_CLASS_EXPORT(LBMKernelETD3Q27CCLBWithSpongeLayer)
BOOST_CLASS_EXPORT(D3Q27EsoTwist3DSplittedVector)
BOOST_CLASS_EXPORT(BCProcessor)
BOOST_CLASS_EXPORT(D3Q27ETBCProcessor)
BOOST_CLASS_EXPORT(D3Q27ETForThinWallBCProcessor)
BOOST_CLASS_EXPORT(DataSet3D)
BOOST_CLASS_EXPORT(Interactor3D)
BOOST_CLASS_EXPORT(D3Q27Interactor)
BOOST_CLASS_EXPORT(CoProcessor)
BOOST_CLASS_EXPORT(MacroscopicQuantitiesCoProcessor)
BOOST_CLASS_EXPORT(ShearStressCoProcessor)
BOOST_CLASS_EXPORT(AverageValuesCoProcessor)
BOOST_CLASS_EXPORT(WbWriterVtkXmlASCII)
BOOST_CLASS_EXPORT(WbWriterVtkXmlBinary)

BOOST_CLASS_EXPORT(BoundaryCondition)
BOOST_CLASS_EXPORT(VelocityBoundaryCondition)
BOOST_CLASS_EXPORT(NonEqDensityBoundaryCondition)
BOOST_CLASS_EXPORT(EqDensityBoundaryCondition)
BOOST_CLASS_EXPORT(NoSlipBoundaryCondition)
BOOST_CLASS_EXPORT(HighViscosityNoSlipBoundaryCondition)
BOOST_CLASS_EXPORT(SlipBoundaryCondition)
BOOST_CLASS_EXPORT(NonReflectingDensityBoundaryCondition)
BOOST_CLASS_EXPORT(NonReflectingVelocityBoundaryCondition)
BOOST_CLASS_EXPORT(ThinWallNoSlipBoundaryCondition)

BOOST_CLASS_EXPORT(D3Q27BoundaryConditionAdapter)
BOOST_CLASS_EXPORT(D3Q27DensityBCAdapter)
BOOST_CLASS_EXPORT(D3Q27NoSlipBCAdapter)
BOOST_CLASS_EXPORT(D3Q27SlipBCAdapter)
BOOST_CLASS_EXPORT(D3Q27VelocityBCAdapter)
#endif
