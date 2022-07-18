#include "LiggghtsCouplingCoProcessor.h"
#include "GbSphere3D.h"
#include "MPICommunicator.h"
#include "CoProcessor.h"
#include "LiggghtsCouplingWrapper.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "DistributionArray3D.h"
#include "DataSet3D.h"
#include "IBcumulantK17LBMKernel.h"

LiggghtsCouplingCoProcessor::LiggghtsCouplingCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                         SPtr<Communicator> comm, LiggghtsCouplingWrapper &wrapper,
                                                         int demSteps)
    : CoProcessor(grid, s), comm(comm), wrapper(wrapper), demSteps(demSteps)
{
    //gridRank     = comm->getProcessID();
    //minInitLevel = this->grid->getCoarsestInitializedLevel();
    //maxInitLevel = this->grid->getFinestInitializedLevel();

    //blockVector.resize(maxInitLevel + 1);

    //for (int level = minInitLevel; level <= maxInitLevel; level++) {
    //    grid->getBlocks(level, gridRank, true, blockVector[level]);
    //}
}

LiggghtsCouplingCoProcessor::~LiggghtsCouplingCoProcessor()
{
}

void LiggghtsCouplingCoProcessor::process(double actualTimeStep)
{ 
    setSpheresOnLattice();
    wrapper.run(demSteps);
    std::cout << "step: " << actualTimeStep << "\n";
}

void LiggghtsCouplingCoProcessor::setSpheresOnLattice()
{
    std::vector<int> excludeType;

    int nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;

    for (int iS = 0; iS < nPart; iS++) 
    {
        int type = (int)wrapper.lmp->atom->type[iS];
        bool excludeFlag(false);
        for (int iT = 0; iT < excludeType.size(); iT++) {
            //std::cout << iS << " " << type << " " << excludeType[iT] << std::endl;
            if (type == excludeType[iT]) {
                excludeFlag = true;
                break;
            }
        }

        if (excludeFlag)
            continue;

        double x[3], v[3], omega[3];
        double r;
        int id = wrapper.lmp->atom->tag[iS];

        for (int i = 0; i < 3; i++) 
        {
            x[i] = wrapper.lmp->atom->x[iS][i];
            v[i]     = wrapper.lmp->atom->v[iS][i];
            omega[i] = wrapper.lmp->atom->omega[iS][i];
        }
        
        r = wrapper.lmp->atom->radius[iS];

        std::cout << "x[0] = " << x[0] << ", x[1] = " << x[1] << ", x[2] = " << x[2] << std::endl;
        std::cout << "v[0] = " << v[0] << ", v[1] = " << v[1] << ", v[2] = " << v[2] << std::endl;
        std::cout << "omega[0] = " << omega[0] << ", omega[1] = " << omega[1] << ", omega[2] = " << omega[2] << std::endl;
        std::cout << "r = " << r << std::endl;
        
        setSingleSphere3D(x, v, omega, r, id);
    }
}

void LiggghtsCouplingCoProcessor::getForcesFromLattice() {}

void LiggghtsCouplingCoProcessor::setSingleSphere3D(double *x, double *v, double *omega, /* double *com,*/ double r,
                                                    int id /*, bool initVelFlag*/)
{
    int level = 0;
    //UbTupleInt3 bi = grid->getBlockIndexes(x[0], x[1], x[2], level);
    //SPtr<Block3D> block = grid->getBlock(val<1>(bi), val<2>(bi), val<3>(bi), level);
    
    std::vector<SPtr<Block3D>> blocks;
    grid->getBlocksByCuboid(level, x[0] - r, x[1] - r, x[2] - r, x[0] + r, x[1] + r, x[2] + r, blocks);
    



    for (SPtr<Block3D> block : blocks) {
        if (block) {
            SPtr<ILBMKernel> kernel = block->getKernel();
            // SPtr<BCArray3D> bcArray                 = kernel->getBCProcessor()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

            CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr particleData =
                dynamicPointerCast<IBcumulantK17LBMKernel>(kernel)->getParticleData();

            if (!particleData)
                continue;

            int minX1 = 1;
            int minX2 = 1;
            int minX3 = 1;

            int maxX1 = (int)(distributions->getNX1()) - 1;
            int maxX2 = (int)(distributions->getNX2()) - 1;
            int maxX3 = (int)(distributions->getNX3()) - 1;

            for (int ix3 = minX3; ix3 < maxX3; ix3++) {
                for (int ix2 = minX2; ix2 < maxX2; ix2++) {
                    for (int ix1 = minX1; ix1 < maxX1; ix1++) {

                        Vector3D nX = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                        
                        double const dx = nX[0] - x[0];
                        double const dy = nX[1] - x[1];
                        double const dz = nX[2] - x[2];

                        double const sf = calcSolidFraction(dx, dy, dz, r);

                        double const sf_old = (*particleData)(ix1,ix2,ix3)->solidFraction;
                        int const id_old = (int)(*particleData)(ix1,ix2,ix3)->partId;

                        int const decFlag = (sf > SOLFRAC_MIN) + 2 * (sf_old > SOLFRAC_MIN);

                        switch (decFlag) {
                            case 0: // sf == 0 && sf_old == 0
                                setToZero(*(*particleData)(ix1, ix2, ix3).get());
                                break; // do nothing
                            case 1:    // sf > 0 && sf_old == 0
                                setValues(*(*particleData)(ix1, ix2, ix3).get(), sf, dx, dy, dz, omega, id);
                                break;
                            case 2:               // sf == 0 && sf_old > 0
                                if (id_old == id) // then particle has left this cell
                                    setToZero(*(*particleData)(ix1, ix2, ix3).get());
                                break; // else do nothing
                            case 3:    // sf > 0 && sf_old > 0
                                if (sf > sf_old || id_old == id)
                                    setValues(*(*particleData)(ix1, ix2, ix3).get(), sf, dx, dy, dz, omega, id);
                                break; // else do nothing
                        }
                        // if desired, initialize interior of sphere with sphere velocity
                       // if (initVelFlag && sf > SOLFRAC_MAX)
                       //     cell.defineVelocity(particleData->uPart);
                    }
                }
            }
        }
    }

}

double LiggghtsCouplingCoProcessor::calcSolidFraction(double const dx_, double const dy_, double const dz_,
                                                      double const r_)
{
    static int const slicesPerDim = 5;
    static double const sliceWidth       = 1. / ((double)slicesPerDim);
    static double const fraction         = 1. / ((double)(slicesPerDim * slicesPerDim * slicesPerDim));

    // should be sqrt(3.)/2.
    // add a little to avoid roundoff errors
    static const double sqrt3half = (double)sqrt(3.1) / 2.;

    double const dist = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;

    double const r_p = r_ + sqrt3half;
    if (dist > r_p * r_p)
        return 0;

    double const r_m = r_ - sqrt3half;
    if (dist < r_m * r_m)
        return 1;

    double const r_sq = r_ * r_;
    double dx_sq[slicesPerDim], dy_sq[slicesPerDim], dz_sq[slicesPerDim];

    // pre-calculate d[xyz]_sq for efficiency
    for (int i = 0; i < slicesPerDim; i++) {
        double const delta = -0.5 + ((double)i + 0.5) * sliceWidth;
        double const dx    = dx_ + delta;
        dx_sq[i]      = dx * dx;
        double const dy    = dy_ + delta;
        dy_sq[i]      = dy * dy;
        double const dz    = dz_ + delta;
        dz_sq[i]      = dz * dz;
    }

    unsigned int n(0);
    for (int i = 0; i < slicesPerDim; i++) {
        for (int j = 0; j < slicesPerDim; j++) {
            for (int k = 0; k < slicesPerDim; k++) {
                n += (dx_sq[i] + dy_sq[j] + dz_sq[k] < r_sq);
            }
        }
    }

    return fraction * ((double)n);
}

  void LiggghtsCouplingCoProcessor::setValues(IBdynamicsParticleData &p, double const sf, double const dx,
                                                 double const dy, double const dz, double* omega, int id)
{
    //p.uPart.from_cArray(v);
    if (omega != 0) {
        p.uPart[0] += omega[1] * dz - omega[2] * dy;
        p.uPart[1] += -omega[0] * dz + omega[2] * dx;
        p.uPart[2] += omega[0] * dy - omega[1] * dx;
    }
    p.solidFraction = sf;
    p.partId        = id;
}


void LiggghtsCouplingCoProcessor::setToZero(IBdynamicsParticleData &p)
{
    p.uPart[0]      = 0;
    p.uPart[1]      = 0;
    p.uPart[2]      = 0;
    p.solidFraction = 0;
    p.partId        = 0;
}