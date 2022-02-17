#include "WritePeBlocksCoProcessor.h"

#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "Block3D.h"
#include <mpi/Communicator.h>
#include "D3Q27System.h"
#include "Grid3D.h"
#include "UbScheduler.h"

WritePeBlocksCoProcessor::WritePeBlocksCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                   WbWriter *const writer, std::shared_ptr<vf::mpi::Communicator> comm,
                                                   SPtr<walberla::blockforest::BlockForest> forest)
    : CoProcessor(grid, s), path(path), writer(writer), comm(comm), forest(forest)
{
}

WritePeBlocksCoProcessor::~WritePeBlocksCoProcessor() {}

void WritePeBlocksCoProcessor::process(double step)
{
    if (scheduler->isDue(step))
        collectData(step);
}

void WritePeBlocksCoProcessor::collectData(double step)
{
    if (comm->getProcessID() == comm->getRoot()) {
        int istep = int(step);
        std::vector<std::string> filenames;
        std::vector<UbTupleFloat3> nodes;
        std::vector<UbTupleInt8> cells;
        std::vector<std::string> celldatanames;

        celldatanames.push_back("ID");
        celldatanames.push_back("rank");

        walberla::uint_t rank = 0;

        std::vector<std::vector<double>> celldata(celldatanames.size());

        int nr = 0;

        for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
            walberla::AABB aabb = blockIt->getAABB();

            nodes.push_back(makeUbTuple((float)aabb.xMin(), (float)aabb.yMin(), (float)aabb.zMin()));
            nodes.push_back(makeUbTuple((float)aabb.xMax(), (float)aabb.yMin(), (float)aabb.zMin()));
            nodes.push_back(makeUbTuple((float)aabb.xMax(), (float)aabb.yMax(), (float)aabb.zMin()));
            nodes.push_back(makeUbTuple((float)aabb.xMin(), (float)aabb.yMax(), (float)aabb.zMin()));
            nodes.push_back(makeUbTuple((float)aabb.xMin(), (float)aabb.yMin(), (float)aabb.zMax()));
            nodes.push_back(makeUbTuple((float)aabb.xMax(), (float)aabb.yMin(), (float)aabb.zMax()));
            nodes.push_back(makeUbTuple((float)aabb.xMax(), (float)aabb.yMax(), (float)aabb.zMax()));
            nodes.push_back(makeUbTuple((float)aabb.xMin(), (float)aabb.yMax(), (float)aabb.zMax()));
            cells.push_back(makeUbTuple(nr, nr + 1, nr + 2, nr + 3, nr + 4, nr + 5, nr + 6, nr + 7));
            nr += 8;

            // data
            celldata[0].push_back((double)blockIt->getId().getID());
            forest->getProcessRank(rank, blockIt->getId());
            celldata[1].push_back((double)rank);
        }

        filenames.push_back(writer->writeOctsWithCellData(
            path + "/peBlocks/peBlocks_" + UbSystem::toString(grid->getRank()) + "_" + UbSystem::toString(istep), nodes,
            cells, celldatanames, celldata));

        if (istep == CoProcessor::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(path + "/peBlocks/peBlocks_collection", filenames,
                                                                istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path + "/peBlocks/peBlocks_collection", filenames,
                                                                     istep, false);
        }

        UBLOG(logINFO, "WritePeBlocksCoProcessor step: " << istep);
    }
}