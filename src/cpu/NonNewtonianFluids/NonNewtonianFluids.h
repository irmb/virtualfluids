#pragma once

#include "NonNewtonianFluids/BoundaryConditions/ThixotropyDensityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNonReflectingOutflowBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityWithDensityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyBinghamModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyHerschelBulkleyModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyPowellEyringModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyBinghamModelVelocityBCStrategy.h"

#include "NonNewtonianFluids/SimulationObservers/CalculateTorqueSimulationObserver.h"
#include "NonNewtonianFluids/SimulationObservers/WriteThixotropyQuantitiesSimulationObserver.h"

#include "NonNewtonianFluids/LBM/ThixotropyLBMKernel.h"
#include "NonNewtonianFluids/LBM/ThixotropyExpLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyBinghamModelLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyHerschelBulkleyModelLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyInterpolator.h"
#include "NonNewtonianFluids/LBM/Rheology.h"
#include "NonNewtonianFluids/LBM/RheologyK17LBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyPowellEyringModelLBMKernel.h"
