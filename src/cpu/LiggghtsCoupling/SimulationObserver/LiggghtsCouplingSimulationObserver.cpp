#include "LiggghtsCouplingSimulationObserver.h"
#include "GbSphere3D.h"
#include "mpi/MPICommunicator.h"
#include "SimulationObserver.h"
#include "LiggghtsCoupling/3rdParty/LiggghtsCouplingWrapper.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "DistributionArray3D.h"
#include "DataSet3D.h"
#include "LiggghtsCoupling/LBM/IBcumulantK17LBMKernel.h"
#include "LiggghtsCoupling/LBM/IBsharpInterfaceLBMKernel.h"
#include "LBMUnitConverter.h"
#include "fix_lb_coupling_onetoone.h"

LiggghtsCouplingSimulationObserver::LiggghtsCouplingSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                         SPtr<vf::parallel::Communicator> comm,
                                                         LiggghtsCouplingWrapper &wrapper, int demSteps,
                                                         SPtr<LBMUnitConverter> units)
    : SimulationObserver(grid, s), comm(comm), wrapper(wrapper), demSteps(demSteps), units(units)
{

}

LiggghtsCouplingSimulationObserver::~LiggghtsCouplingSimulationObserver()
{
}

void LiggghtsCouplingSimulationObserver::update(double actualTimeStep)
{ 
    //if (comm->getProcessID() == 0)
    //    std::cout << "LiggghtsCouplingSimulationObserver step: " << actualTimeStep << "\n";
    
    //comm->barrier();

    getForcesFromLattice();

    //comm->barrier();
    
    wrapper.run(demSteps);

    //comm->barrier();
    
    setSpheresOnLattice();

    //comm->barrier();
}

void LiggghtsCouplingSimulationObserver::setSpheresOnLattice()
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

        double x[3] = { 0, 0, 0 }, v[3] = { 0, 0, 0 }, omega[3] = { 0, 0, 0 };
        double r;
        int id = wrapper.lmp->atom->tag[iS];

        for (int i = 0; i < 3; i++) 
        {
            x[i]     = wrapper.lmp->atom->x[iS][i]; // * units->getFactorLentghWToLb(); // - 0.5; ????
            v[i]     = wrapper.lmp->atom->v[iS][i] * units->getFactorVelocityWToLb();
            omega[i] = wrapper.lmp->atom->omega[iS][i] / units->getFactorTimeWToLb();
        }
        
        r = wrapper.lmp->atom->radius[iS]; // * units->getFactorLentghWToLb();

        //std::cout << "x[0] = " << x[0] << ", x[1] = " << x[1] << ", x[2] = " << x[2] << std::endl;
        //std::cout << "v[0] = " << v[0] << ", v[1] = " << v[1] << ", v[2] = " << v[2] << std::endl;
        //std::cout << "omega[0] = " << omega[0] << ", omega[1] = " << omega[1] << ", omega[2] = " << omega[2] << std::endl;
        //std::cout << "r = " << r << std::endl;
        
        setSingleSphere3D(x, v, omega, r, id);
    }
}

void LiggghtsCouplingSimulationObserver::setSingleSphere3D(double *x, double *v, double *omega, /* double *com,*/ double r,
                                                    int id /*, bool initVelFlag*/)
{
    int level = 0;
    //UbTupleInt3 bi = grid->getBlockIndexes(x[0], x[1], x[2], level);
    //SPtr<Block3D> block = grid->getBlock(val<1>(bi), val<2>(bi), val<3>(bi), level);
    
    std::vector<SPtr<Block3D>> blocks;
    grid->getBlocksByCuboid(level, x[0] - r, x[1] - r, x[2] - r, x[0] + r, x[1] + r, x[2] + r, blocks);

    //DEBUG
    ///////////////////////
    if (blocks.size() == 2) 
        int test = 0;

    ///////////////////////

    for (SPtr<Block3D> block : blocks) {
        if (block) {
            SPtr<ILBMKernel> kernel = block->getKernel();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

            CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr particleData =
                dynamicPointerCast<LiggghtsCouplingLBMKernel>(kernel)->getParticleData();

            if (!particleData)
                continue;

            int minX1b = 1;
            int minX2b = 1;
            int minX3b = 1;

            int maxX1b = (int)(distributions->getNX1()) - 2;
            int maxX2b = (int)(distributions->getNX2()) - 2;
            int maxX3b = (int)(distributions->getNX3()) - 2;

            real deltax = grid->getDeltaX(block);

            UbTupleInt3 nodesMin = grid->getNodeIndexes(block, x[0] - r - deltax, x[1] - r - deltax, x[2] - r - deltax);
            UbTupleInt3 nodesMax = grid->getNodeIndexes(block, x[0] + r + deltax, x[1] + r + deltax, x[2] + r + deltax);

            int minX1 = (val<1>(nodesMin) < minX1b) ? minX1b : val<1>(nodesMin);
            int minX2 = (val<2>(nodesMin) < minX2b) ? minX2b : val<2>(nodesMin);
            int minX3 = (val<3>(nodesMin) < minX3b) ? minX3b : val<3>(nodesMin);

            int maxX1 = (val<1>(nodesMax) > maxX1b) ? maxX1b : val<1>(nodesMax);
            int maxX2 = (val<2>(nodesMax) > maxX2b) ? maxX2b : val<2>(nodesMax);
            int maxX3 = (val<3>(nodesMax) > maxX3b) ? maxX3b : val<3>(nodesMax);

            //int minX1 =  minX1b;
            //int minX2 =  minX2b;
            //int minX3 =  minX3b;

            //int maxX1 =  maxX1b;
            //int maxX2 =  maxX2b;
            //int maxX3 =  maxX3b;


            for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
                for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                    for (int ix1 = minX1; ix1 <= maxX1; ix1++) {

                        //UbTupleInt3 blockNX = grid->getBlockNX();

                        //double const dx = val<1>(blockNX) * block->getX1() + ix1 - x[0];
                        //double const dy = val<2>(blockNX) * block->getX2() + ix2 - x[1];
                        //double const dz = val<3>(blockNX) * block->getX3() + ix3 - x[2];

                        Vector3D worldCoordinates = grid->getNodeCoordinates(block, ix1, ix2, ix3);

                        double const dx = (worldCoordinates[0] - x[0]) * units->getFactorLentghWToLb();
                        double const dy = (worldCoordinates[1] - x[1]) * units->getFactorLentghWToLb();
                        double const dz = (worldCoordinates[2] - x[2]) * units->getFactorLentghWToLb();

                        double const sf = calcSolidFraction(dx, dy, dz, r * units->getFactorLentghWToLb());

                        double const sf_old = (*particleData)(ix1,ix2,ix3)->solidFraction;
                        int const id_old = (int)(*particleData)(ix1,ix2,ix3)->partId;

                        int const decFlag = (sf > SOLFRAC_MIN) + 2 * (sf_old > SOLFRAC_MIN);

                        switch (decFlag) {
                            case 0: // sf == 0 && sf_old == 0
                                setToZero(*(*particleData)(ix1, ix2, ix3).get());
                                break; // do nothing
                            case 1:    // sf > 0 && sf_old == 0
                                setValues(*(*particleData)(ix1, ix2, ix3).get(), id, sf, v, dx, dy, dz, omega);
                                break;
                            case 2:               // sf == 0 && sf_old > 0
                                if (id_old == id) // then particle has left this cell
                                    setToZero(*(*particleData)(ix1, ix2, ix3).get());
                                break; // else do nothing
                            case 3:    // sf > 0 && sf_old > 0
                                if (sf > sf_old || id_old == id)
                                    setValues(*(*particleData)(ix1, ix2, ix3).get(), id, sf, v, dx, dy, dz, omega);
                                break; // else do nothing
                        }
                        // if desired, initialize interior of sphere with sphere velocity
                       // if (initVelFlag && sf > SOLFRAC_MAX)
                       //     cell.defineVelocity(particleData->uPart);

                        //if (sf > 0) {
                        //    std::cout << "sf = " << sf << std::endl;
                        //    std::cout << "ix1 = " << ix1 << ", ix2 = " << ix2 << ", ix3 = " << ix3 << std::endl;
                        //}
                    }
                }
            }
        }
    }

}

double LiggghtsCouplingSimulationObserver::calcSolidFraction(double const dx_, double const dy_, double const dz_,
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
    double dx_sq[slicesPerDim] = { 0, 0, 0, 0, 0 }, dy_sq[slicesPerDim] = { 0, 0, 0, 0, 0 }, dz_sq[slicesPerDim] = { 0, 0, 0, 0, 0 };

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

  void LiggghtsCouplingSimulationObserver::setValues(IBdynamicsParticleData &p, int const id, double const sf, double const *v, double const dx, double const dy, double const dz, double const *omega)
{
    p.uPart[0] = v[0];
    p.uPart[1] = v[1];
    p.uPart[2] = v[2];

    if (omega != 0) {
        p.uPart[0] += omega[1] * dz - omega[2] * dy;
        p.uPart[1] += -omega[0] * dz + omega[2] * dx;
        p.uPart[2] += omega[0] * dy - omega[1] * dx;
    }
    p.solidFraction = sf;
    p.partId        = id;
}


void LiggghtsCouplingSimulationObserver::setToZero(IBdynamicsParticleData &p)
{
    p.uPart[0]      = 0;
    p.uPart[1]      = 0;
    p.uPart[2]      = 0;
    p.solidFraction = 0;
    p.partId        = 0;
}

void LiggghtsCouplingSimulationObserver::getForcesFromLattice()
{
    static std::vector<double> force, torque;
    static typename ParticleData::ParticleDataArrayVector x_lb;

    int const nPart   = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    int const n_force = nPart * 3;

    if (nPart == 0)
        return; // no particles - no work

    if (nPart > (int)x_lb.size()) {
        for (int iPart = 0; iPart < (int)x_lb.size(); iPart++) {
            x_lb[iPart][0] = wrapper.lmp->atom->x[iPart][0];
            x_lb[iPart][1] = wrapper.lmp->atom->x[iPart][1];
            x_lb[iPart][2] = wrapper.lmp->atom->x[iPart][2];
        }
        for (int iPart = (int)x_lb.size(); iPart < nPart; iPart++) {
            std::array<double, 3> ar = {wrapper.lmp->atom->x[iPart][0],
                                        wrapper.lmp->atom->x[iPart][1],
                                        wrapper.lmp->atom->x[iPart][2]};
            x_lb.push_back(ar);
        }
            

    } else {
        for (int iPart = 0; iPart < nPart; iPart++) {
            x_lb[iPart][0] = wrapper.lmp->atom->x[iPart][0];
            x_lb[iPart][1] = wrapper.lmp->atom->x[iPart][1];
            x_lb[iPart][2] = wrapper.lmp->atom->x[iPart][2];
        }
    }

    if (n_force > (int)force.size()) {
        for (int i = 0; i < (int)force.size(); i++) {
            force[i]  = 0;
            torque[i] = 0;
        }
        for (int i = (int)force.size(); i < n_force; i++) {
            force.push_back(0.);
            torque.push_back(0.);
        }
    } else {
        for (int i = 0; i < n_force; i++) {
            force[i]  = 0;
            torque[i] = 0;
        }
    }

    SumForceTorque3D(x_lb, &force.front(), &torque.front());

    LAMMPS_NS::FixLbCouplingOnetoone *couplingFix =
        dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone *>(wrapper.lmp->modify->find_fix_style("couple/lb/onetoone", 0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    for (int iPart = 0; iPart < nPart; iPart++)
        for (int j = 0; j < 3; j++) {
            f_liggghts[iPart][j] = 0;
            t_liggghts[iPart][j] = 0;
        }

    for (int iPart = 0; iPart < nPart; iPart++) {
        int tag          = wrapper.lmp->atom->tag[iPart];
        int liggghts_ind = wrapper.lmp->atom->map(tag);

        for (int j = 0; j < 3; j++) {
            f_liggghts[liggghts_ind][j] += force[3 * iPart + j] * units->getFactorForceLbToW();
            t_liggghts[liggghts_ind][j] += torque[3 * iPart + j] * units->getFactorTorqueLbToW();
        }
    }
    couplingFix->comm_force_torque();
}

void LiggghtsCouplingSimulationObserver::SumForceTorque3D(ParticleData::ParticleDataArrayVector &x, double *force, double *torque)
{
    int nx = grid->getNX1(), ny = grid->getNX2(), nz = grid->getNX3();

    std::vector < SPtr < Block3D > > blocks;
    int level = 0;
    grid->getBlocks(level, grid->getRank(), true, blocks);

        
    for (SPtr<Block3D> block : blocks) {
        if (block) {
            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

            CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr particleData =
                dynamicPointerCast<LiggghtsCouplingLBMKernel>(kernel)->getParticleData();

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

                        // LIGGGHTS indices start at 1
                        int const id = (*particleData)(ix1, ix2, ix3)->partId;
                        if (id < 1)
                            continue; // no particle here

                        int const ind = wrapper.lmp->atom->map(id);

                        if (ind < 0) continue; // no particle here

                        Vector3D worldCoordinates = grid->getNodeCoordinates(block, ix1, ix2, ix3);

                        double dx = (worldCoordinates[0] - x[ind][0]) * units->getFactorLentghWToLb();
                        double dy = (worldCoordinates[1] - x[ind][1]) * units->getFactorLentghWToLb();
                        double dz = (worldCoordinates[2] - x[ind][2]) * units->getFactorLentghWToLb();

                        // minimum image convention, needed if
                        // (1) PBC are used and
                        // (2) both ends of PBC lie on the same processor
                        if ((int)dx > nx / 2)
                            dx -= nx;
                        else if ((int)dx < -nx / 2)
                            dx += nx;
                        if ((int)dy > ny / 2)
                            dy -= ny;
                        else if ((int)dy < -ny / 2)
                            dy += ny;
                        if ((int)dz > nz / 2)
                            dz -= nz;
                        else if ((int)dz < -nz / 2)
                            dz += nz;

                        double const forceX = (*particleData)(ix1, ix2, ix3)->hydrodynamicForce[0];
                        double const forceY = (*particleData)(ix1, ix2, ix3)->hydrodynamicForce[1];
                        double const forceZ = (*particleData)(ix1, ix2, ix3)->hydrodynamicForce[2];

                        double const torqueX = dy * forceZ - dz * forceY;
                        double const torqueY = -dx * forceZ + dz * forceX;
                        double const torqueZ = dx * forceY - dy * forceX;

                        addForce(ind, 0, forceX, force);
                        addForce(ind, 1, forceY, force);
                        addForce(ind, 2, forceZ, force);

                        addTorque(ind, 0, torqueX, torque);
                        addTorque(ind, 1, torqueY, torque);
                        addTorque(ind, 2, torqueZ, torque);
                    }
                }
            }
        }
    }
 }

void LiggghtsCouplingSimulationObserver::addForce(int const partId, int const coord, double const value, double *force)
{
    force[3 * partId + coord] += value;
}

void LiggghtsCouplingSimulationObserver::addTorque(int const partId, int const coord, double const value, double *torque)
{
    torque[3 * partId + coord] += value;
}