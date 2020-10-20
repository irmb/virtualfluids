#include "CalculateForcesCoProcessor.h"
#include "BCProcessor.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Communicator.h"
#include "D3Q27Interactor.h"
#include "DataSet3D.h"
#include "DistributionArray3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

CalculateForcesCoProcessor::CalculateForcesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                       SPtr<Communicator> comm, double v, double a)
    : CoProcessor(grid, s), path(path), comm(comm), v(v), a(a), forceX1global(0), forceX2global(0), forceX3global(0)
{
    if (comm->getProcessID() == comm->getRoot()) {
        std::ofstream ostr;
        std::string fname = path;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string file_path = UbSystem::getPathFromString(fname);
            if (file_path.size() > 0) {
                UbSystem::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr.width(12);
        ostr << "step"
             << "\t";
        ostr.width(12);
        ostr << "Cx"
             << "\t";
        ostr.width(12);
        ostr << "Cy"
             << "\t";
        ostr.width(12);
        ostr << "Cz"
             << "\t";
        ostr.width(12);
        ostr << "Fx"
             << "\t";
        ostr.width(12);
        ostr << "Fy"
             << "\t";
        ostr.width(12);
        ostr << "Fz" << std::endl;
        ostr.close();
    }
}
//////////////////////////////////////////////////////////////////////////
CalculateForcesCoProcessor::~CalculateForcesCoProcessor() = default;
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::process(double step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "D3Q27ForcesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::collectData(double step)
{
    calculateForces();

    if (comm->getProcessID() == comm->getRoot()) {
        int istep = static_cast<int>(step);
        std::ofstream ostr;
        std::string fname = path;
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

        calculateCoefficients();

        ostr.width(12);
        ostr.setf(std::ios::fixed);
        ostr << istep << "\t";
        write(&ostr, C1, (char *)"\t");
        write(&ostr, C2, (char *)"\t");
        write(&ostr, C3, (char *)"\t");
        write(&ostr, forceX1global, (char *)"\t");
        write(&ostr, forceX2global, (char *)"\t");
        write(&ostr, forceX3global, (char *)"\t");
        ostr << std::endl;
        ostr.close();
    }
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::calculateForces()
{
    forceX1global = 0.0;
    forceX2global = 0.0;
    forceX3global = 0.0;

    for (SPtr<D3Q27Interactor> interactor : interactors) {
        for (BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
            double forceX1 = 0.0;
            double forceX2 = 0.0;
            double forceX3 = 0.0;

            SPtr<Block3D> block                             = t.first;
            std::set<std::vector<int>> &transNodeIndicesSet = t.second;

            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCProcessor()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            distributions->swap();

            int ghostLayerWidth = kernel->getGhostLayerWidth();
            int minX1           = ghostLayerWidth;
            int maxX1           = (int)bcArray->getNX1() - 1 - ghostLayerWidth;
            int minX2           = ghostLayerWidth;
            int maxX2           = (int)bcArray->getNX2() - 1 - ghostLayerWidth;
            int minX3           = ghostLayerWidth;
            int maxX3           = (int)bcArray->getNX3() - 1 - ghostLayerWidth;

            for (std::vector<int> node : transNodeIndicesSet) {
                int x1 = node[0];
                int x2 = node[1];
                int x3 = node[2];

                // without ghost nodes
                if (x1 < minX1 || x1 > maxX1 || x2 < minX2 || x2 > maxX2 || x3 < minX3 || x3 > maxX3)
                    continue;

                if (bcArray->isFluid(
                        x1, x2,
                        x3)) // es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
                {
                    SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
                    UbTupleDouble3 forceVec     = getForces(x1, x2, x3, distributions, bc);
                    forceX1 += val<1>(forceVec);
                    forceX2 += val<2>(forceVec);
                    forceX3 += val<3>(forceVec);
                }
            }
            // if we have got discretization with more level
            // deltaX is LBM deltaX and equal LBM deltaT
            double deltaX = LBMSystem::getDeltaT(block->getLevel()); // grid->getDeltaT(block);
            double deltaXquadrat = deltaX * deltaX;
            forceX1 *= deltaXquadrat;
            forceX2 *= deltaXquadrat;
            forceX3 *= deltaXquadrat;

            distributions->swap();

            forceX1global += forceX1;
            forceX2global += forceX2;
            forceX3global += forceX3;
        }
    }
    std::vector<double> values;
    std::vector<double> rvalues;
    values.push_back(forceX1global);
    values.push_back(forceX2global);
    values.push_back(forceX3global);

    rvalues = comm->gather(values);
    if (comm->getProcessID() == comm->getRoot()) {
        forceX1global = 0.0;
        forceX2global = 0.0;
        forceX3global = 0.0;

        for (int i = 0; i < (int)rvalues.size(); i += 3) {
            forceX1global += rvalues[i];
            forceX2global += rvalues[i + 1];
            forceX3global += rvalues[i + 2];
        }
    }
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 CalculateForcesCoProcessor::getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions,
                                                     SPtr<BoundaryConditions> bc)
{
    UbTupleDouble3 force(0.0, 0.0, 0.0);

    if (bc) {
        // references to tuple "force"
        double &forceX1 = val<1>(force);
        double &forceX2 = val<2>(force);
        double &forceX3 = val<3>(force);
        double f, fnbr;

        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            if (bc->hasNoSlipBoundaryFlag(fdir)) {
                const int invDir = D3Q27System::INVDIR[fdir];
                f = dynamicPointerCast<EsoTwist3D>(distributions)->getDistributionInvForDirection(x1, x2, x3, invDir);
                fnbr =
                    dynamicPointerCast<EsoTwist3D>(distributions)
                        ->getDistributionInvForDirection(x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir],
                                                         x3 + D3Q27System::DX3[invDir], fdir);

                forceX1 += (f + fnbr) * D3Q27System::DX1[invDir];
                forceX2 += (f + fnbr) * D3Q27System::DX2[invDir];
                forceX3 += (f + fnbr) * D3Q27System::DX3[invDir];
            }
        }
    }

    return force;
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::calculateCoefficients()
{
    double F1 = forceX1global;
    double F2 = forceX2global;
    double F3 = forceX3global;

    // return 2*F/(rho*v*v*a);
    C1 = 2.0 * F1 / (v * v * a);
    C2 = 2.0 * F2 / (v * v * a);
    C3 = 2.0 * F3 / (v * v * a);
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::addInteractor(SPtr<D3Q27Interactor> interactor) { interactors.push_back(interactor); }
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::write(std::ofstream *fileObject, double value, char *separator)
{
    (*fileObject).width(12);
    //(*fileObject).precision(2);
    (*fileObject).setf(std::ios::fixed);
    (*fileObject) << value;
    (*fileObject) << separator;
}
