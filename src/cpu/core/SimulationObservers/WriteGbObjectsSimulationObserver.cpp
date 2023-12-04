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
//! \file WriteGbObjectsSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#include "WriteGbObjectsSimulationObserver.h"
#include <parallel/Communicator.h>
#include "GbObject3D.h"
#include "UbScheduler.h"
#include "WbWriterVtkXmlASCII.h"
#include "WbWriterVtkXmlBinary.h"
#include <vector>

WriteGbObjectsSimulationObserver::WriteGbObjectsSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                     WbWriter *const writer, std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), writer(writer), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
WriteGbObjectsSimulationObserver::~WriteGbObjectsSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsSimulationObserver::addGbObject(SPtr<GbObject3D> object) { objects.push_back(object); }
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsSimulationObserver::collectData(real step)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt3> triangles;

    int numObjcts = 0;

    for (SPtr<GbObject3D> object : objects) {
        object->addSurfaceTriangleSet(nodes, triangles);
        numObjcts++;
    }

    int istep = static_cast<int>(step);

    std::string pfilePath, partPath, subfolder, cfilePath;

    subfolder = "gob" + UbSystem::toString(istep);
    pfilePath = path + "/gob/" + subfolder;
    cfilePath = path + "/gob/gob_collection";
    partPath  = pfilePath + "/gob" + UbSystem::toString(comm->getProcessID()) + "_" + UbSystem::toString(istep);

    std::string partName = writer->writeTriangles(partPath, nodes, triangles);

    size_t found      = partName.find_last_of("/");
    std::string piece = partName.substr(found + 1);
    piece             = subfolder + "/" + piece;

    std::vector<std::string> datanames;
    std::vector<std::string> cellDataNames;
    // std::vector<std::string> pieces = comm->gather(piece);
    std::vector<std::string> pieces;
    pieces.push_back(piece);
    if (comm->isRoot()) {
        std::string pname =
            WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        found = pname.find_last_of("/");
        piece = pname.substr(found + 1);

        std::vector<std::string> filenames;
        filenames.push_back(piece);
        if (step == SimulationObserver::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "WriteGbObjectsSimulationObserver number of objects: " << numObjcts);
        UBLOG(logINFO, "WriteGbObjectsSimulationObserver step: " << istep);
    }
}
