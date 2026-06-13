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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>
#include <functional>

#include "gpu/GridGenerator/global.h"

#include <basics/geometry3d/GbSpatialData3D.h>

#include "gpu/GridGenerator/grid/NodeValues.h"

namespace vf::gpu {

class Grid;

class Side;
enum class SideType;

class TransientBCInputFileReader;

namespace grid_generator
{
class BoundaryCondition
{
public:
    virtual ~BoundaryCondition() = default;

    std::vector<uint> indices;
    SPtr<Side> side;
    std::vector<std::vector<real>> qs;

    std::vector<uint> patches;

    virtual char getType() const = 0;

    bool isSide(SideType side) const;

    real getQ(uint index, uint dir) { return this->qs[index][dir]; }

    void getCoords( SPtr<Grid> grid, std::vector<real>& x, std::vector<real>& y, std::vector<real>& z);
    virtual void setAdditionalIndices(const SPtr<Grid>& /*grid*/, uint /*index*/) {};

};

}

//////////////////////////////////////////////////////////////////////////

class PressureBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<PressureBoundaryCondition> make(real rho)
    {
        return SPtr<PressureBoundaryCondition>(new PressureBoundaryCondition(rho));
    }

    // matrix indices!!!
    std::vector<uint> neighborIndices;

    real rho;
protected:
    PressureBoundaryCondition(real rho) : rho(rho) { }

public:
    char getType() const override
    {
        return vf::gpu::BC_PRESSURE;
    }

    real getRho() const 
    {
        return this->rho;
    }
    void setAdditionalIndices(const SPtr<Grid> &grid, uint index) override;
};

//////////////////////////////////////////////////////////////////////////

class SlipBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<SlipBoundaryCondition> make(real normalX, real normalY, real normalZ)
    {
        return SPtr<SlipBoundaryCondition>(new SlipBoundaryCondition(normalX, normalY, normalZ));
    }

    real normalX, normalY, normalZ;
    std::vector<real> normalXList, normalYList, normalZList;
protected:
    SlipBoundaryCondition(real normalX, real normalY, real normalZ) : normalX(normalX), normalY(normalY), normalZ(normalZ) { }

public:
    char getType() const override
    {
        return vf::gpu::BC_SLIP;
    }

    void fillSlipNormalLists()
    {   
        for (uint index : this->indices) {
            (void)index;
            this->normalXList.push_back(normalX);
            this->normalYList.push_back(normalY);
            this->normalZList.push_back(normalZ);
        }
    }

    real getNormalx() const { return this->normalX; }
    real getNormaly() const { return this->normalY; }
    real getNormalz() const { return this->normalZ; }

    real getNormalx(uint index) { return this->normalXList[index]; }
    real getNormaly(uint index) { return this->normalYList[index]; }
    real getNormalz(uint index) { return this->normalZList[index]; }
};

//////////////////////////////////////////////////////////////////////////

class StressBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<StressBoundaryCondition> make(real normalX, real normalY, real normalZ, uint samplingOffset,
                                              real vonKarmanConstant, real roughnessLength,
                                              std::shared_ptr<GbSpatialData3D<real>> roughnessMap)
    {
        return SPtr<StressBoundaryCondition>(new StressBoundaryCondition(
            normalX, normalY, normalZ, samplingOffset, vonKarmanConstant, roughnessLength, std::move(roughnessMap)));
    }

    const real normalX, normalY, normalZ;
    const uint samplingOffset;
    const real vonKarmanConstant, roughnessLength;
    std::shared_ptr<GbSpatialData3D<real>> roughnessMap;
    std::vector<real> normalXList, normalYList, normalZList;
    std::vector<real> samplingDistanceList;
    std::vector<real> roughnessLengthList;
    std::vector<uint> samplingIndices;

protected:
    StressBoundaryCondition(real normalX, real normalY, real normalZ, uint samplingOffset, real vonKarmanConstant,
                            real roughnessLength, std::shared_ptr<GbSpatialData3D<real>> roughnessMap)
        : normalX(normalX), normalY(normalY), normalZ(normalZ), samplingOffset(samplingOffset),
          vonKarmanConstant(vonKarmanConstant), roughnessLength(roughnessLength), roughnessMap(std::move(roughnessMap))
    {
    }

public:
    char getType() const override
    {
        return vf::gpu::BC_STRESS;
    }
    
    void fillLists()
    {
        std::fill_n(std::back_inserter(this->normalXList), this->indices.size(), normalX);
        std::fill_n(std::back_inserter(this->normalYList), this->indices.size(), normalY);
        std::fill_n(std::back_inserter(this->normalZList), this->indices.size(), normalZ);
        if(roughnessMap == nullptr)
            std::fill_n(std::back_inserter(this->roughnessLengthList), this->indices.size(), roughnessLength);
    }

    real getNormalx() const { return this->normalX; }
    real getNormaly() const { return this->normalY; }
    real getNormalz() const { return this->normalZ; }

    real getNormalx(uint index) const { return this->normalXList[index]; }
    real getNormaly(uint index) const { return this->normalYList[index]; }
    real getNormalz(uint index) const { return this->normalZList[index]; }

    uint getSamplingOffset() const { return this->samplingOffset; }
    uint getSamplingIndex(uint index) const { return this->samplingIndices[index]; }
    real getRoughnessLength() const { return this->roughnessLength; }
    real getRoughnessLength(uint index) const { return this->roughnessLengthList[index]; }
    real getSamplingDistance(uint index) const { return this->samplingDistanceList[index]; }
    real getVonKarmanConstant() const {return vonKarmanConstant;}

    void setAdditionalIndices(const SPtr<Grid>& grid, uint index) override;

};

class SurfaceLayerBoundaryCondition : public StressBoundaryCondition
{
public:
    static SPtr<SurfaceLayerBoundaryCondition> make(real normalX, real normalY, real normalZ, uint samplingOffset,
                                                    real vonKarmanConstant, real roughnessLength,
                                                    real roughnessLengthTemperature, real surfaceHeatFlux,
                                                    real surfaceTemperature, real heatingRate,
                                                    std::shared_ptr<GbSpatialData3D<real>> roughnessMap)
    {
        return SPtr<SurfaceLayerBoundaryCondition>(new SurfaceLayerBoundaryCondition(
            normalX, normalY, normalZ, samplingOffset, vonKarmanConstant, roughnessLength, roughnessLengthTemperature,
            surfaceHeatFlux, surfaceTemperature, heatingRate, std::move(roughnessMap)));
    }

    real roughnessLengthTemperature, surfaceHeatFlux, surfaceTemperature, heatingRate;
    std::vector<real> roughnessLengthTemperatureList, surfaceHeatFluxList, surfaceTemperatureList, heatingRateList;

protected:
    SurfaceLayerBoundaryCondition(real normalX, real normalY, real normalZ, uint samplingOffset, real vonKarmanConstant,
                                  real roughnessLength, real roughnessLengthTemperature, real surfaceHeatFlux,
                                  real surfaceTemperature, real heatingRate,
                                  std::shared_ptr<GbSpatialData3D<real>> roughnessMap)
        : StressBoundaryCondition(normalX, normalY, normalZ, samplingOffset, vonKarmanConstant, roughnessLength,
                                  std::move(roughnessMap)),
          roughnessLengthTemperature(roughnessLengthTemperature), surfaceHeatFlux(surfaceHeatFlux),
          surfaceTemperature(surfaceTemperature), heatingRate(heatingRate)
    {
    }

public:
    void fillLists()
    {
        StressBoundaryCondition::fillLists();
        if(roughnessMap)
            roughnessLengthTemperatureList = roughnessLengthList;
        else
            std::fill_n(std::back_inserter(roughnessLengthTemperatureList), indices.size(), roughnessLengthTemperature);
        std::fill_n(std::back_inserter(surfaceHeatFluxList), indices.size(), surfaceHeatFlux);
        std::fill_n(std::back_inserter(surfaceTemperatureList), indices.size(), surfaceTemperature);
        std::fill_n(std::back_inserter(heatingRateList), indices.size(), heatingRate);
        std::fill_n(std::back_inserter(roughnessLengthTemperatureList), indices.size(), roughnessLengthTemperature);
    }

    real getRoughnessLengthTemperature() const
    {
        return this->roughnessLengthTemperature;
    }
    real getRoughnessLengthTemperature(uint index)
    {
        return this->roughnessLengthTemperatureList[index];
    }
    real getSurfaceHeatFlux() const
    {
        return this->surfaceHeatFlux;
    }
    real getSurfaceHeatFlux(uint index)
    {
        return this->surfaceHeatFluxList[index];
    }
    real getSurfaceTemperature() const
    {
        return this->surfaceTemperature;
    }
    real getSurfaceTemperature(uint index)
    {
        return this->surfaceTemperatureList[index];
    }
    real getHeatingRate() const
    {
        return this->heatingRate;
    }
    real getHeatingRate(uint index)
    {
        return this->heatingRateList[index];
    }
};

//////////////////////////////////////////////////////////////////////////

class VelocityBoundaryCondition : public grid_generator ::BoundaryCondition
{
public:
    static SPtr<VelocityBoundaryCondition> make(real vx, real vy, real vz)
    {
        return SPtr<VelocityBoundaryCondition>(new VelocityBoundaryCondition(vx, vy, vz));
    }

    real vx, vy, vz;
    std::vector<real> vxList, vyList, vzList;
protected:
    VelocityBoundaryCondition(real vx, real vy, real vz) : vx(vx), vy(vy), vz(vz) { }

public:
    char getType() const override
    {
        return vf::gpu::BC_VELOCITY;
    }

    void fillVelocityLists()
    {
        for( uint index : this->indices ) {
            (void) index;
            this->vxList.push_back(vx);
            this->vyList.push_back(vy);
            this->vzList.push_back(vz);
        }
    }

    real getVx() const { return this->vx; }
    real getVy() const { return this->vy; }
    real getVz() const { return this->vz; }

    real getVx(uint index) { return this->vxList[index]; }
    real getVy(uint index) { return this->vyList[index]; }
    real getVz(uint index) { return this->vzList[index]; }


    void setVelocityProfile( SPtr<Grid> grid, std::function<void(real,real,real,real&,real&,real&)> velocityProfile );
};

//////////////////////////////////////////////////////////////////////////


class GeometryBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<GeometryBoundaryCondition> make()
    {
        return SPtr<GeometryBoundaryCondition>(new GeometryBoundaryCondition());
    }

    real vx, vy, vz;

    real normalX, normalY, normalZ;
    std::vector<real> normalXList, normalYList, normalZList;

    std::vector<real> vxList, vyList, vzList;

private:
    GeometryBoundaryCondition(real vx = 0.0, real vy = 0.0, real vz = 0.0, real normalX = 0.0, real normalY = 0.0, real normalZ = 0.0) :
        vx(vx), vy(vy), vz(vz), normalX(normalX), normalY(normalY), normalZ(normalZ) { }

public:
    char getType() const override
    {
        return vf::gpu::BC_SOLID;
    }

    void setVelocityForPatch( uint patch, real vx, real vy, real vz ){
        for( uint index = 0; index < this->indices.size(); index++ ){
            if( this->patches[index] == patch ){
                this->vxList[index] = vx;
                this->vyList[index] = vy;
                this->vzList[index] = vz;
            }
        }
    }

    void setTangentialVelocityForPatch( SPtr<Grid> grid, uint patch,
                                                  real p1x, real p1y, real p1z, 
                                                  real p2x, real p2y, real p2z, 
                                                  real v, real r );

    void fillVelocityLists()
    {
        for( uint index : this->indices ) {
            (void) index;
            this->vxList.push_back(vx);
            this->vyList.push_back(vy);
            this->vzList.push_back(vz);
        }
    }

    real getVx() const { return this->vx; }
    real getVy() const { return this->vy; }
    real getVz() const { return this->vz; }

    real getVx(uint index) { return this->vxList[index]; }
    real getVy(uint index) { return this->vyList[index]; }
    real getVz(uint index) { return this->vzList[index]; }


    void fillSlipNormalLists()
    {   
        for (uint index : this->indices) {
            (void)index;
            this->normalXList.push_back(normalX);
            this->normalYList.push_back(normalY);
            this->normalZList.push_back(normalZ);
        }
    }

    real getNormalx() const { return this->normalX; }
    real getNormaly() const { return this->normalY; }
    real getNormalz() const { return this->normalZ; }

    real getNormalx(uint index) { return this->normalXList[index]; }
    real getNormaly(uint index) { return this->normalYList[index]; }
    real getNormalz(uint index) { return this->normalZList[index]; }
};

class PrecursorBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<PrecursorBoundaryCondition> make(SPtr<TransientBCInputFileReader> reader, int timeStepsBetweenReads, real velocityX, real velocityY, real velocityZ)
    {
        return SPtr<PrecursorBoundaryCondition>(new PrecursorBoundaryCondition(reader, timeStepsBetweenReads, velocityX, velocityY, velocityZ));
    }

    SPtr<TransientBCInputFileReader> getReader(){ return reader; }
    real getVelocityX() const { return velocityX; }
    real getVelocityY() const { return velocityY; }
    real getVelocityZ() const { return velocityZ; }

private:
    PrecursorBoundaryCondition(SPtr<TransientBCInputFileReader> _reader, uint _timeStepsBetweenReads, real vx, real vy, real vz) : timeStepsBetweenReads(_timeStepsBetweenReads), velocityX(vx), velocityY(vy), velocityZ(vz), reader(_reader) { };
    char getType() const override
    {
        return vf::gpu::BC_VELOCITY;
    }
public:
    uint timeStepsBetweenReads; //!> read data every nth timestep

private:
    real velocityX = 0.0;
    real velocityY = 0.0;
    real velocityZ = 0.0;
    SPtr<TransientBCInputFileReader> reader;
};

class ADNoFluxBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<ADNoFluxBoundaryCondition> make()
    {
        return SPtr<ADNoFluxBoundaryCondition>(new ADNoFluxBoundaryCondition());
    }

protected:
    ADNoFluxBoundaryCondition() = default;

public:
    char getType() const override
    {
        return vf::gpu::BC_AD;
    }
};


class ADFluxBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<ADFluxBoundaryCondition> make(real normalX, real normalY, real normalZ, real gradient)
    {
        return SPtr<ADFluxBoundaryCondition>(new ADFluxBoundaryCondition(normalX, normalY, normalZ, gradient));
    }

    real normalX, normalY, normalZ, gradient;
    std::vector<real> normalXList, normalYList, normalZList, gradientList;

protected:
    ADFluxBoundaryCondition(real normalX, real normalY, real normalZ, real gradient) : normalX(normalX), normalY(normalY), normalZ(normalZ), gradient(gradient)
    {
    }

public:
    char getType() const override
    {
        return vf::gpu::BC_AD;
    }

    void fillBoundaryValueLists();

    real getNormalX(uint index) {return this->normalXList[index];}
    real getNormalY(uint index) {return this->normalYList[index];}
    real getNormalZ(uint index) {return this->normalZList[index];}
    real getGradient(uint index) {return this->gradientList[index];}
};


class ADDirichletBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<ADDirichletBoundaryCondition> make(real BCvalue, real vx, real vy, real vz)
    {
        return SPtr<ADDirichletBoundaryCondition>(new ADDirichletBoundaryCondition(BCvalue, vx, vy, vz));
    }

protected:
    ADDirichletBoundaryCondition(real BCValue, real vx, real vy, real vz) : BCvalue(BCValue), vx(vx), vy(vy), vz(vz) 
    {

    }

public:
    char getType() const override
    {
        return vf::gpu::BC_AD;
    }

    void fillBoundaryValueLists()
    {
        std::fill_n(std::back_inserter(this->vxList), this->indices.size(), vx);
        std::fill_n(std::back_inserter(this->vyList), this->indices.size(), vy);
        std::fill_n(std::back_inserter(this->vzList), this->indices.size(), vz);
        std::fill_n(std::back_inserter(this->BCvalueList), this->indices.size(), BCvalue);
    }

    real getBCvalue() const { return this->BCvalue; }
    real getVx() const { return this->vx; }
    real getVy() const { return this->vy; }
    real getVz() const { return this->vz; }

    real getBCvalue(uint index) { return this->BCvalueList[index]; }
    real getVx(uint index)  { return this->vxList[index]; }
    real getVy(uint index)  { return this->vyList[index]; }
    real getVz(uint index)  { return this->vzList[index]; }

private:
    real BCvalue, vx, vy, vz;
    std::vector<real> BCvalueList, vxList, vyList, vzList;
};

class ADOutflowBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<ADOutflowBoundaryCondition> make()
    {
        return SPtr<ADOutflowBoundaryCondition>(new ADOutflowBoundaryCondition());
    }

    // matrix indices of the neighbor node in the outflow direction (kN)
    std::vector<uint> neighborIndices;

protected:
    ADOutflowBoundaryCondition() = default;

public:
    char getType() const override
    {
        return vf::gpu::BC_AD;
    }

    void setAdditionalIndices(const SPtr<Grid>& grid, uint index) override;
};


class ADNeumannBoundaryCondition : public grid_generator::BoundaryCondition
{
public:
    static SPtr<ADNeumannBoundaryCondition> make(real BCvalue, real vx, real vy, real vz)
    {
        return SPtr<ADNeumannBoundaryCondition>(new ADNeumannBoundaryCondition(BCvalue, vx, vy, vz));
    }

    real gradient, vx, vy, vz;
    std::vector<real> gradientList, vxList, vyList, vzList;

protected:
    ADNeumannBoundaryCondition(real gradient, real vx, real vy, real vz) : gradient(gradient), vx(vx), vy(vy), vz(vz) 
    {

    }

public:
    char getType() const override
    {
        return vf::gpu::BC_AD;
    }

    void fillBoundaryValueLists();

    real getBCgradient() const { return this->gradient; }
    real getVx() const { return this->vx; }
    real getVy() const { return this->vy; }
    real getVz() const { return this->vz; }

    real getBCgradient(uint index) { return this->gradientList[index]; }
    real getVx(uint index)  { return this->vxList[index]; }
    real getVy(uint index)  { return this->vyList[index]; }
    real getVz(uint index)  { return this->vzList[index]; }
};

}

#endif
//! \}
