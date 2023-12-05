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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GPU_INTERFACE_H
#define GPU_INTERFACE_H

#include "LBM/LB.h"

#include <cuda.h>
#include <cuda_runtime.h>

struct LBMSimulationParameter;
class Parameter;


void QADDev7(unsigned int numberOfThreads,
                        real* DD, 
                        real* DD7,
                        real* temp,
                        real diffusivity,
                        int* k_Q, 
                        real* QQ,
                        unsigned int numberOfBCnodes, 
                        real om1, 
                        unsigned int* neighborX,
                        unsigned int* neighborY,
                        unsigned int* neighborZ,
                        unsigned long long numberOfLBnodes, 
                        bool isEvenTimestep);

//////////////////////////////////////////////////////////////////////////
//! \brief defines the behavior of a slip-AD boundary condition
void ADSlipVelDevComp(
    uint numberOfThreads,
    real * normalX,
    real * normalY,
    real * normalZ,
    real * distributions,
    real * distributionsAD,
    int* QindexArray,
    real * Qarrays,
    uint numberOfBCnodes,
    real omegaDiffusivity,
    uint * neighborX,
    uint * neighborY,
    uint * neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);
    
void QADDirichletDev27( unsigned int numberOfThreads,
                                   real* DD, 
                                   real* DD27,
                                   real* temp,
                                   real diffusivity,
                                   int* k_Q, 
                                   real* QQ,
                                   unsigned int numberOfBCnodes, 
                                   real om1, 
                                   unsigned int* neighborX,
                                   unsigned int* neighborY,
                                   unsigned int* neighborZ,
                                   unsigned long long numberOfLBnodes, 
                                   bool isEvenTimestep);

void QADBBDev27(  unsigned int numberOfThreads,
                             real* DD, 
                             real* DD27,
                             real* temp,
                             real diffusivity,
                             int* k_Q, 
                             real* QQ,
                             unsigned int numberOfBCnodes, 
                             real om1, 
                             unsigned int* neighborX,
                             unsigned int* neighborY,
                             unsigned int* neighborZ,
                             unsigned long long numberOfLBnodes, 
                             bool isEvenTimestep);

void QADVelDev7(unsigned int numberOfThreads,
                           real* DD, 
                           real* DD7,
                           real* temp,
                           real* velo,
                           real diffusivity,
                           int* k_Q, 
                           real* QQ,
                           unsigned int numberOfBCnodes, 
                           real om1, 
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           unsigned long long numberOfLBnodes, 
                           bool isEvenTimestep);


void QADVelDev27(  unsigned int numberOfThreads,
                              real* DD, 
                              real* DD27,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q, 
                              real* QQ,
                              unsigned int numberOfBCnodes, 
                              real om1, 
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned long long numberOfLBnodes, 
                              bool isEvenTimestep);

void QADPressDev7( unsigned int numberOfThreads,
                              real* DD, 
                              real* DD7,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q, 
                              real* QQ,
                              unsigned int numberOfBCnodes, 
                              real om1, 
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned long long numberOfLBnodes, 
                              bool isEvenTimestep);

void QADPressDev27(unsigned int numberOfThreads,
                              real* DD, 
                              real* DD27,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q, 
                              real* QQ,
                              unsigned int numberOfBCnodes, 
                              real om1, 
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned long long numberOfLBnodes, 
                              bool isEvenTimestep);

void QADPressNEQNeighborDev27(
                                            unsigned int numberOfThreads,
                                            real* DD,
                                            real* DD27,
                                            int* k_Q,
                                            int* k_N,
                                            int numberOfBCnodes,
                                            unsigned int* neighborX,
                                            unsigned int* neighborY,
                                            unsigned int* neighborZ,
                                            unsigned long long numberOfLBnodes,
                                            bool isEvenTimestep
                                        );

void QNoSlipADincompDev7(unsigned int numberOfThreads,
                                    real* DD, 
                                    real* DD7,
                                    real* temp,
                                    real diffusivity,
                                    int* k_Q, 
                                    real* QQ,
                                    unsigned int numberOfBCnodes, 
                                    real om1, 
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned long long numberOfLBnodes, 
                                    bool isEvenTimestep);

void QNoSlipADincompDev27(unsigned int numberOfThreads,
                                     real* DD, 
                                     real* DD27,
                                     real* temp,
                                     real diffusivity,
                                     int* k_Q, 
                                     real* QQ,
                                     unsigned int numberOfBCnodes, 
                                     real om1, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned long long numberOfLBnodes, 
                                     bool isEvenTimestep);

void QADVeloIncompDev7( unsigned int numberOfThreads,
                                   real* DD, 
                                   real* DD7,
                                   real* temp,
                                   real* velo,
                                   real diffusivity,
                                   int* k_Q, 
                                   real* QQ,
                                   unsigned int numberOfBCnodes, 
                                   real om1, 
                                   unsigned int* neighborX,
                                   unsigned int* neighborY,
                                   unsigned int* neighborZ,
                                   unsigned long long numberOfLBnodes, 
                                   bool isEvenTimestep);


void QADVeloIncompDev27( unsigned int numberOfThreads,
                                    real* DD, 
                                    real* DD27,
                                    real* temp,
                                    real* velo,
                                    real diffusivity,
                                    int* k_Q, 
                                    real* QQ,
                                    unsigned int numberOfBCnodes, 
                                    real om1, 
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned long long numberOfLBnodes, 
                                    bool isEvenTimestep);

void QADPressIncompDev7(  unsigned int numberOfThreads,
                                     real* DD, 
                                     real* DD7,
                                     real* temp,
                                     real* velo,
                                     real diffusivity,
                                     int* k_Q, 
                                     real* QQ,
                                     unsigned int numberOfBCnodes, 
                                     real om1, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned long long numberOfLBnodes, 
                                     bool isEvenTimestep);

void QADPressIncompDev27(  unsigned int numberOfThreads,
                                      real* DD, 
                                      real* DD27,
                                      real* temp,
                                      real* velo,
                                      real diffusivity,
                                      int* k_Q, 
                                      real* QQ,
                                      unsigned int numberOfBCnodes, 
                                      real om1, 
                                      unsigned int* neighborX,
                                      unsigned int* neighborY,
                                      unsigned int* neighborZ,
                                      unsigned long long numberOfLBnodes, 
                                      bool isEvenTimestep);

#endif
