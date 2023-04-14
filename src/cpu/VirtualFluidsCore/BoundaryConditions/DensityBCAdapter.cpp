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
//! \file DensityBCAdapter.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "DensityBCAdapter.h"
#include "basics/utilities/UbInfinity.h"
#include "basics/utilities/UbLogger.h"

using namespace std;
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const real &dens, const real &startTime, const real &endTime)
{
    this->densBCs.emplace_back(dens, startTime, endTime);
    this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const BCFunction &densBC)
{
    this->densBCs.push_back(densBC);
    this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const std::vector<BCFunction> &densBCs)
{
    this->densBCs = densBCs;
    this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const mu::Parser &function, const real &startTime, const real &endTime)
{
    this->densBCs.emplace_back(function, startTime, endTime);
    this->init();
}
/*==========================================================*/
void DensityBCAdapter::init()
{
    this->timeStep = 0.0;

    this->x1 = 0.0;
    this->x2 = 0.0;
    this->x3 = 0.0;

    this->tmpDensityFunction = NULL;

    try // initilialization and validation of functions
    {
        for (size_t pos = 0; pos < densBCs.size(); ++pos) {
            if (!(UbMath::equal(BCFunction::INFCONST, densBCs[pos].getEndTime()) &&
                  UbMath::greaterEqual(this->timeStep, densBCs[pos].getStartTime()))) {
                this->setTimeDependent();
            }

            densBCs[pos].getFunction().DefineVar("t", &this->timeStep);
            densBCs[pos].getFunction().DefineVar("x1", &this->x1);
            densBCs[pos].getFunction().DefineVar("x2", &this->x2);
            densBCs[pos].getFunction().DefineVar("x3", &this->x3);

            densBCs[pos].getFunction().Eval(); //<-- validation
        }
    } catch (mu::Parser::exception_type &e) {
        stringstream error;
        error << "mu::parser exception occurs, message(" << e.GetMsg() << "), formula("
              << e.GetExpr() + "), token(" + e.GetToken() << ")"
              << ", pos(" << e.GetPos() << "), error code(" << e.GetCode();
        throw UbException(error.str());
    } catch (...) {
        throw UbException(UB_EXARGS, "unknown exception");
    }
}
/*==========================================================*/
void DensityBCAdapter::init(const D3Q27Interactor *const & /*interactor*/, const real &time)
{
    this->timeStep           = time;
    this->tmpDensityFunction = NULL;
    real maxEndtime        = -Ub::inf;

    // aktuelle Densityfunction bestimmen
    for (size_t pos = 0; pos < densBCs.size(); ++pos) {
        if (UbMath::equal(densBCs[pos].getEndTime(), BCFunction::INFTIMEDEPENDENT))
            maxEndtime = Ub::inf;
        maxEndtime = UbMath::max(maxEndtime, densBCs[pos].getStartTime(),
                                 densBCs[pos].getEndTime()); // startTime abfragen, da  INFCONST=-10

        if (UbMath::greaterEqual(this->timeStep, densBCs[pos].getStartTime())) {
            if (UbMath::lessEqual(this->timeStep, densBCs[pos].getEndTime()) ||
                UbMath::equal(densBCs[pos].getEndTime(), (real)BCFunction::INFCONST) ||
                UbMath::equal(densBCs[pos].getEndTime(), (real)BCFunction::INFTIMEDEPENDENT)) {
                tmpDensityFunction = &densBCs[pos].getFunction();
                break;
            }
        }
    }

    // wenn funktionen zweitlich konstant sind und bis t=unendlich gelten
    // kann man zeitabhaengigkeit deaktivieren
    if (UbMath::greaterEqual(time, maxEndtime))
        this->unsetTimeDependent();

    UBLOG(logDEBUG4, "D3Q27DensityBCAdapter::init(time="
                         << time << ") "
                         << ", rho= \"" << (tmpDensityFunction ? tmpDensityFunction->GetExpr() : "-")
                         << "\", timedependant=" << (this->isTimeDependent() ? "true" : "false"));
}
/*==========================================================*/
void DensityBCAdapter::update(const D3Q27Interactor *const &interactor, const real &time)
{
    this->init(interactor, time);
}
/*==========================================================*/
void DensityBCAdapter::adaptBCForDirection(const D3Q27Interactor & /*interactor*/, SPtr<BoundaryConditions> bc,
                                           const real & /*worldX1*/, const real & /*worldX2*/,
                                           const real & /*worldX3*/, const real &q, const int &fdirection,
                                           const real & /*time*/)
{
    bc->setDensityBoundaryFlag(D3Q27System::INVDIR[fdirection], secondaryBcOption);
    bc->setQ((real)q, fdirection);
}
/*==========================================================*/
void DensityBCAdapter::adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                               const real &worldX2, const real &worldX3, const real &time)
{
    this->setNodeDensity(interactor, bc, worldX1, worldX2, worldX3, time);
    bc->setBcAlgorithmType(algorithmType);
}
/*==========================================================*/
void DensityBCAdapter::setNodeDensity(const D3Q27Interactor & /*interactor*/, SPtr<BoundaryConditions> bc,
                                      const real &worldX1, const real &worldX2, const real &worldX3,
                                      const real &timestep)
{
    // Geschwindigkeiten setzen
    try {
        // PunktKoordinaten bestimmen
        this->x1       = worldX1;
        this->x2       = worldX2;
        this->x3       = worldX3;
        this->timeStep = timestep;

        if (tmpDensityFunction)
            bc->setBoundaryDensity((real)tmpDensityFunction->Eval());
    } catch (mu::Parser::exception_type &e) {
        stringstream error;
        error << "mu::parser exception occurs, message(" << e.GetMsg() << "), formula("
              << e.GetExpr() + "), token(" + e.GetToken() << ")"
              << ", pos(" << e.GetPos() << "), error code(" << e.GetCode();
        throw UbException(error.str());
    } catch (...) {
        throw UbException(UB_EXARGS, "unknown exception");
    }
}
/*==========================================================*/
real DensityBCAdapter::getDensity(const real &x1, const real &x2, const real &x3, const real &timeStep)
{
    this->x1       = x1;
    this->x2       = x2;
    this->x3       = x3;
    this->timeStep = timeStep;

    if (!tmpDensityFunction)
        return 0.0;

    return tmpDensityFunction->Eval();
}
/*==========================================================*/
string DensityBCAdapter::toString()
{
    stringstream info;
    info << "D3Q27DensityBCAdapter:\n";
    info << " #dens-functions = " << (int)densBCs.size() << endl;
    info << " protected variables: x1, x2, x3, t" << endl;

    for (size_t i = 0; i < densBCs.size(); ++i) {
        info << "\n   dens-function nr." << i << ":" << endl;
        info << densBCs[i].toString() << endl;
    }
    return info.str();
}
