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
//! \file WriteMultiphaseQuantitiesCoProcessor.h
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef WriteThixotropyQuantitiesCoProcessor_H
#define WriteThixotropyQuantitiesCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include <mpi/Communicator.h>
#include "WbWriter.h"

class WriteThixotropyQuantitiesCoProcessor : public  CoProcessor
{
public:
	WriteThixotropyQuantitiesCoProcessor();
	WriteThixotropyQuantitiesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<LBMUnitConverter> conv, std::shared_ptr<vf::mpi::Communicator> comm);
	~WriteThixotropyQuantitiesCoProcessor() = default;

   void process(real step) override;

protected:
   void collectData(real step);
   void addDataMQ(SPtr<Block3D> block);
   void clearData();

private:
   void init();
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleUInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<real> > data; 
   std::string path;
   WbWriter* writer;
   SPtr<LBMUnitConverter> conv;
//   bool bcInformation;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   std::shared_ptr<vf::mpi::Communicator> comm;
//	double ConcentrationSum;
};
#endif
