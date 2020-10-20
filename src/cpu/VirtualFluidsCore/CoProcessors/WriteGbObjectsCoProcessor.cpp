#include "WriteGbObjectsCoProcessor.h"
#include "Communicator.h"
#include "GbObject3D.h"
#include "UbScheduler.h"
#include "WbWriterVtkXmlASCII.h"
#include "WbWriterVtkXmlBinary.h"
#include <vector>

WriteGbObjectsCoProcessor::WriteGbObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                     WbWriter *const writer, SPtr<Communicator> comm)
    : CoProcessor(grid, s), path(path), writer(writer), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
WriteGbObjectsCoProcessor::~WriteGbObjectsCoProcessor() = default;
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsCoProcessor::process(double step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsCoProcessor::addGbObject(SPtr<GbObject3D> object) { objects.push_back(object); }
//////////////////////////////////////////////////////////////////////////
void WriteGbObjectsCoProcessor::collectData(double step)
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
        if (step == CoProcessor::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "WriteGbObjectsCoProcessor number of objects: " << numObjcts);
        UBLOG(logINFO, "WriteGbObjectsCoProcessor step: " << istep);
    }
}
