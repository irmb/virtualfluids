#include "LineTimeSeriesCoProcessor.h"
#include "BCProcessor.h"
#include "WbWriterVtkXmlASCII.h"

#include "Block3D.h"
#include <mpi/Communicator.h>
#include "CompressibleCumulantLBMKernel.h"
#include "CoordinateTransformation3D.h"
#include "DataSet3D.h"
#include "GbLine3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

LineTimeSeriesCoProcessor::LineTimeSeriesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                     SPtr<GbLine3D> line, int level, std::shared_ptr<vf::mpi::Communicator> comm)
    : CoProcessor(grid, s), path(path), length(0), ix1(0), ix2(0), ix3(0), level(level), line(line)
{
    root  = comm->isRoot();
    fname = path;

    mpi_comm  = *((MPI_Comm *)comm->getNativeCommunicator());
    numOfProc = comm->getNumberOfProcesses();
    gridRank  = comm->getProcessID();

    double dx = CoProcessor::grid->getDeltaX(level);

    SPtr<CoordinateTransformation3D> trafo = grid->getCoordinateTransformator();
    double orgX1                           = trafo->getX1CoordinateOffset();
    double orgX2                           = trafo->getX2CoordinateOffset();
    double orgX3                           = trafo->getX3CoordinateOffset();

    int x1min = (int)((line->getX1Minimum() - orgX1) / dx);
    int x1max = (int)((line->getX1Maximum() - orgX1) / dx);
    int x2min = (int)((line->getX2Minimum() - orgX2) / dx);
    int x2max = (int)((line->getX2Maximum() - orgX2) / dx);
    int x3min = (int)((line->getX3Minimum() - orgX3) / dx);
    int x3max = (int)((line->getX3Maximum() - orgX3) / dx);

    UbTupleInt3 blockNx = grid->getBlockNX();

    if (x1min != x1max) {
        dir     = X1;
        blocknx = val<1>(blockNx);
        length  = x1max;
    } else if (x2min != x2max) {
        dir     = X2;
        blocknx = val<2>(blockNx);
        length  = x2max;
    } else if (x3min != x3max) {
        dir     = X3;
        blocknx = val<3>(blockNx);
        length  = x3max;
    }
    blockix1 = x1min / val<1>(blockNx);
    blockix2 = x2min / val<2>(blockNx);
    blockix3 = x3min / val<3>(blockNx);

    ix1 = x1min % val<1>(blockNx) + 1;
    ix2 = x2min % val<2>(blockNx) + 1;
    ix3 = x3min % val<3>(blockNx) + 1;
}
//////////////////////////////////////////////////////////////////////////
void LineTimeSeriesCoProcessor::process(double step)
{
    if (scheduler->isDue(step)) {
        collectData();
    }

    UBLOG(logDEBUG3, "MacroscopicQuantitiesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void LineTimeSeriesCoProcessor::writeLine(const std::string &path)
{
    std::vector<UbTupleFloat3> nodes(2);
    std::vector<UbTupleInt2> lines(1);
    val<1>(nodes[0])            = (float)line->getX1Minimum();
    val<2>(nodes[0])            = (float)line->getX2Minimum();
    val<3>(nodes[0])            = (float)line->getX3Minimum();
    val<1>(nodes[1])            = (float)line->getX1Maximum();
    val<2>(nodes[1])            = (float)line->getX2Maximum();
    val<3>(nodes[1])            = (float)line->getX3Maximum();
    val<1>(lines[0])            = 0;
    val<1>(lines[0])            = 1;
    WbWriterVtkXmlASCII *writer = WbWriterVtkXmlASCII::getInstance();
    writer->writeLines(path, nodes, lines);
}
//////////////////////////////////////////////////////////////////////////
void LineTimeSeriesCoProcessor::collectData()
{
    LBMReal f[27];
    LBMReal vx1, vx2, vx3, rho;
    MPI_Status status;
    std::vector<double> v1(length, 0);
    std::vector<double> v2(length, 0);
    std::vector<double> v3(length, 0);
    std::vector<double> p(length, 0);
    for (int x = 0; x < length; x += blocknx) {
        if (dir == X1) {
            blockix1 = x / blocknx;
        } else if (dir == X2) {
            blockix2 = x / blocknx;
        } else if (dir == X3) {
            blockix3 = x / blocknx;
        }

        SPtr<Block3D> block = CoProcessor::grid->getBlock(blockix1, blockix2, blockix3, level);
        if (block) {
            if (block->getRank() == gridRank) {
                SPtr<ILBMKernel> kernel = block->getKernel();
                calcMacros              = NULL;
                if (kernel->getCompressible()) {
                    calcMacros = &D3Q27System::calcCompMacroscopicValues;
                } else {
                    calcMacros = &D3Q27System::calcIncompMacroscopicValues;
                }
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

                for (int ix = 1; ix <= blocknx; ix++) {
                    if (dir == X1) {
                        ix1 = ix;
                    } else if (dir == X2) {
                        ix2 = ix;
                    } else if (dir == X3) {
                        ix3 = ix;
                    }
                    distributions->getDistribution(f, ix1, ix2, ix3);
                    calcMacros(f, rho, vx1, vx2, vx3);
                    v1[x + (ix - 1)] = vx1;
                    v2[x + (ix - 1)] = vx2;
                    v3[x + (ix - 1)] = vx3;
                    p[x + (ix - 1)]  = rho;
                }
            }
        }
    }

    if (root) {
        for (int i = 1; i < numOfProc; i++) {
            std::vector<double> v1temp(length, 0);
            std::vector<double> v2temp(length, 0);
            std::vector<double> v3temp(length, 0);
            std::vector<double> ptemp(length, 0);
            MPI_Recv(&v1temp[0], length, MPI_DOUBLE, i, 1, mpi_comm, &status);
            MPI_Recv(&v2temp[0], length, MPI_DOUBLE, i, 2, mpi_comm, &status);
            MPI_Recv(&v3temp[0], length, MPI_DOUBLE, i, 3, mpi_comm, &status);
            MPI_Recv(&ptemp[0], length, MPI_DOUBLE, i, 4, mpi_comm, &status);
            for (int j = 0; j < length; j++) {
                v1[j] += v1temp[j];
                v2[j] += v2temp[j];
                v3[j] += v3temp[j];
                p[j] += ptemp[j];
            }
        }

        std::ofstream ostr;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        for (int x = 0; x < length; x++) {
            ostr << v1[x];
            if (x < length - 1) {
                ostr << ";";
            }
        }
        ostr << "\n";
        for (int x = 0; x < length; x++) {
            ostr << v2[x];
            if (x < length - 1) {
                ostr << ";";
            }
        }
        ostr << "\n";
        for (int x = 0; x < length; x++) {
            ostr << v3[x];
            if (x < length - 1) {
                ostr << ";";
            }
        }
        ostr << "\n";
        for (int x = 0; x < length; x++) {
            ostr << p[x];
            if (x < length - 1) {
                ostr << ";";
            }
        }
        ostr << "\n";
        ostr.close();
    } else {
        MPI_Send(&v1[0], length, MPI_DOUBLE, 0, 1, mpi_comm);
        MPI_Send(&v2[0], length, MPI_DOUBLE, 0, 2, mpi_comm);
        MPI_Send(&v3[0], length, MPI_DOUBLE, 0, 3, mpi_comm);
        MPI_Send(&p[0], length, MPI_DOUBLE, 0, 4, mpi_comm);
    }
}
