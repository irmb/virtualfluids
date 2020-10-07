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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include <PointerDefinitions.h>
#include "LBMSystem.h"
#include "ILBMKernel.h"
#include <array>
#include <limits>
#include <muParser.h>

class BCProcessor;
class DataSet3D;
class Block3D;

//! \brief A base class provides basic functionality for LBM kernel 
class LBMKernel : public ILBMKernel, public enableSharedFromThis<LBMKernel>
{
public:
    typedef std::numeric_limits<LBMReal> LBMRealLim;
public:
    LBMKernel();

    virtual SPtr<LBMKernel> clone() = 0;

    void setBCProcessor(SPtr<BCProcessor> bcp) override;
    SPtr<BCProcessor> getBCProcessor() const override;

    void setCollisionFactor(double collFactor) override;
    double getCollisionFactor() const override;

    void setGhostLayerWidth(int witdh);
    int  getGhostLayerWidth() const override;

    void setDataSet(SPtr<DataSet3D> dataSet);
    SPtr<DataSet3D> getDataSet() const override;

    void setForcingX1(LBMReal forcingX1);
    void setForcingX2(LBMReal forcingX2);
    void setForcingX3(LBMReal forcingX3);

    void setForcingX1(const mu::Parser& parser);
    void setForcingX2(const mu::Parser& parser);
    void setForcingX3(const mu::Parser& parser);

    void setForcingX1(const std::string& muParserString);
    void setForcingX2(const std::string& muParserString);
    void setForcingX3(const std::string& muParserString);

    void setIndex(int x1, int x2, int x3);

    LBMReal getDeltaT() const override;
    void setDeltaT(LBMReal dt);

    bool getCompressible() const override;
    void setCompressible(bool val);

    bool getWithForcing() const override;
    void setWithForcing(bool val);

    bool getWithSpongeLayer() const;
    void setWithSpongeLayer(bool val);

    void setSpongeLayer(const mu::Parser& parser);
    void setSpongeLayer(const std::string& muParserString);

    void setBlock(SPtr<Block3D> block);
    SPtr<Block3D> getBlock() const;

    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const override;

    void swapDistributions() override;

    void setNX(std::array<int, 3> nx);
    std::array<int, 3> getNX();

protected:
    SPtr<DataSet3D> dataSet;
    SPtr<BCProcessor> bcProcessor;
    LBMReal collFactor;
    int ghostLayerWidth{1};
    bool compressible{false};

    //forcing 
    bool withForcing{false};
    mu::Parser muForcingX1;
    mu::Parser muForcingX2;
    mu::Parser muForcingX3;
    int ix1, ix2, ix3;
    LBMReal deltaT{1.0};

    //sponge layer
    bool withSpongeLayer{false};
    mu::Parser muSpongeLayer;

    WPtr<Block3D> block;

    std::array<int, 3> nx;

private:
    void checkFunction(mu::Parser fct);
};

#endif
