//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Data Data
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef DataSet3D_h
#define DataSet3D_h

#include <PointerDefinitions.h>

#include "DistributionArray3D.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

using AverageValuesArray3D     = CbArray4D<real, IndexerX4X3X2X1>;
using ShearStressValuesArray3D = CbArray4D<real, IndexerX4X3X2X1>;
using RelaxationFactorArray3D  = CbArray3D<real, IndexerX3X2X1>;
using PhaseFieldArray3D        = CbArray3D<real, IndexerX3X2X1>;
using PressureFieldArray3D     = CbArray3D<real, IndexerX3X2X1>;

//! A class provides an interface for data structures in the kernel.
class DataSet3D
{
public:
    SPtr<DistributionArray3D> getFdistributions() const;
    void setFdistributions(SPtr<DistributionArray3D> distributions);

    SPtr<DistributionArray3D> getHdistributions() const;
    void setHdistributions(SPtr<DistributionArray3D> distributions);

    //SPtr<DistributionArray3D> getH1distributions() const;
    //void setH1distributions(SPtr<DistributionArray3D> distributions);

    SPtr<DistributionArray3D> getH2distributions() const;
    void setH2distributions(SPtr<DistributionArray3D> distributions);

    SPtr<AverageValuesArray3D> getAverageDensity() const;
    void setAverageDensity(SPtr<AverageValuesArray3D> values);

    SPtr<AverageValuesArray3D> getAverageVelocity() const;
    void setAverageVelocity(SPtr<AverageValuesArray3D> values);

    SPtr<AverageValuesArray3D> getAverageFluctuations() const;
    void setAverageFluctuations(SPtr<AverageValuesArray3D> values);

    SPtr<AverageValuesArray3D> getAverageTriplecorrelations() const;
    void setAverageTriplecorrelations(SPtr<AverageValuesArray3D> values);

    SPtr<AverageValuesArray3D> getAverageValues() const;
    void setAverageValues(SPtr<AverageValuesArray3D> values);

    SPtr<ShearStressValuesArray3D> getShearStressValues() const;
    void setShearStressValues(SPtr<ShearStressValuesArray3D> values);

    SPtr<RelaxationFactorArray3D> getRelaxationFactor() const;
    void setRelaxationFactor(SPtr<RelaxationFactorArray3D> values);

    SPtr<PhaseFieldArray3D> getPhaseField() const;
    void setPhaseField(SPtr<PhaseFieldArray3D> values);

    SPtr<PhaseFieldArray3D> getPhaseField2() const;
    void setPhaseField2(SPtr<PhaseFieldArray3D> values);

    SPtr<PressureFieldArray3D> getPressureField() const;
    void setPressureField(SPtr<PressureFieldArray3D> values);

protected:
private:
    SPtr<DistributionArray3D> fdistributions;
    SPtr<DistributionArray3D> hdistributions;
    //SPtr<DistributionArray3D> h1distributions;
    SPtr<DistributionArray3D> h2distributions;
 
    SPtr<AverageValuesArray3D> averageValues;
    SPtr<AverageValuesArray3D> averageDensity;
    SPtr<AverageValuesArray3D> averageVelocity;
    SPtr<AverageValuesArray3D> averageFluktuations;
    SPtr<AverageValuesArray3D> averageTriplecorrelations;
    SPtr<ShearStressValuesArray3D> shearStressValues;

    SPtr<RelaxationFactorArray3D> relaxationFactor;
    
    SPtr<PhaseFieldArray3D> phaseField;
    SPtr<PhaseFieldArray3D> phaseField2;
    SPtr<PressureFieldArray3D> pressureField;
};

inline SPtr<DistributionArray3D> DataSet3D::getFdistributions() const { return fdistributions; }

inline void DataSet3D::setFdistributions(SPtr<DistributionArray3D> distributions) { fdistributions = distributions; }

inline SPtr<DistributionArray3D> DataSet3D::getHdistributions() const { return hdistributions; }

inline void DataSet3D::setHdistributions(SPtr<DistributionArray3D> distributions) { hdistributions = distributions; }

//inline SPtr<DistributionArray3D> DataSet3D::getH1distributions() const { return h1distributions; }
//
//inline void DataSet3D::setH1distributions(SPtr<DistributionArray3D> distributions) { h1distributions = distributions; }

inline SPtr<DistributionArray3D> DataSet3D::getH2distributions() const { return h2distributions; }

inline void DataSet3D::setH2distributions(SPtr<DistributionArray3D> distributions) { h2distributions = distributions; }

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageValues() const { return averageValues; }

inline void DataSet3D::setAverageValues(SPtr<AverageValuesArray3D> values) { averageValues = values; }

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageDensity() const { return averageDensity; }

inline void DataSet3D::setAverageDensity(SPtr<AverageValuesArray3D> values) { averageDensity = values; }

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageVelocity() const { return averageVelocity; }

inline void DataSet3D::setAverageVelocity(SPtr<AverageValuesArray3D> values) { averageVelocity = values; }

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageFluctuations() const { return averageFluktuations; }

inline void DataSet3D::setAverageFluctuations(SPtr<AverageValuesArray3D> values) { averageFluktuations = values; }

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageTriplecorrelations() const { return averageTriplecorrelations; }

inline void DataSet3D::setAverageTriplecorrelations(SPtr<AverageValuesArray3D> values)
{
    averageTriplecorrelations = values;
}

inline SPtr<ShearStressValuesArray3D> DataSet3D::getShearStressValues() const { return shearStressValues; }

inline void DataSet3D::setShearStressValues(SPtr<ShearStressValuesArray3D> values) { shearStressValues = values; }

inline SPtr<RelaxationFactorArray3D> DataSet3D::getRelaxationFactor() const { return relaxationFactor; }

inline void DataSet3D::setRelaxationFactor(SPtr<RelaxationFactorArray3D> values) { relaxationFactor = values; }

inline SPtr<PhaseFieldArray3D> DataSet3D::getPhaseField() const { return phaseField; }

inline void DataSet3D::setPhaseField(SPtr<PhaseFieldArray3D> values) { phaseField = values; }

inline SPtr<PhaseFieldArray3D> DataSet3D::getPhaseField2() const { return phaseField2; }

inline void DataSet3D::setPhaseField2(SPtr<PhaseFieldArray3D> values) { phaseField2 = values; }

inline SPtr<PressureFieldArray3D> DataSet3D::getPressureField() const { return pressureField; }

inline void DataSet3D::setPressureField(SPtr<PressureFieldArray3D> values) { pressureField = values; }

#endif

//! \}
