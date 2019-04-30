#include "GridProvider.h"

#include <Parameter/Parameter.h>
#include "GridReaderFiles/GridReader.h"
#include "GridReaderGenerator/GridGenerator.h"

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <GPU/CudaMemoryManager.h>



void GridProvider::setNumberOfNodes(const int numberOfNodes, const int level) const
{
	para->getParH(level)->size_Mat_SP = numberOfNodes;
	para->getParD(level)->size_Mat_SP = numberOfNodes;
	para->getParH(level)->mem_size_real_SP = sizeof(real) * para->getParH(level)->size_Mat_SP;
	para->getParH(level)->mem_size_int_SP = sizeof(uint) * para->getParH(level)->size_Mat_SP;
	para->getParD(level)->mem_size_real_SP = sizeof(real) * para->getParD(level)->size_Mat_SP;
	para->getParD(level)->mem_size_int_SP = sizeof(uint) * para->getParD(level)->size_Mat_SP;
}

void GridProvider::setInitalNodeValues(const int numberOfNodes, const int level) const
{
	const real PI = 3.141592653589793238462643383279;
	const real gridX = para->getParH(level)->gridNX - 1;
	const real gridY = para->getParH(level)->gridNY - 1;
	const real gridZ = para->getParH(level)->gridNZ - 1;
	////////////////////////////////////////////////////////////////////////////////
	//standard initial conditions
	for (int j = 1; j <= numberOfNodes; j++)
	{
		para->getParH(level)->rho_SP[j] = 0.0;

		para->getParH(level)->vx_SP[j] = para->getVelocity();// 0.0;
		para->getParH(level)->vy_SP[j] = 0.0;
		para->getParH(level)->vz_SP[j] = 0.0;

		if (para->getCalcMedian()) {
			para->getParH(level)->vx_SP_Med[j] = 0.0f;
			para->getParH(level)->vy_SP_Med[j] = 0.0f;
			para->getParH(level)->vz_SP_Med[j] = 0.0f;
			para->getParH(level)->rho_SP_Med[j] = 0.0f;
			para->getParH(level)->press_SP_Med[j] = 0.0f;
		}
		if (para->getUseWale()) {
			para->getParH(level)->turbViscosity[j] = 0.0f;
			//Debug
			para->getParH(level)->gSij[j] = 0.0f;
			para->getParH(level)->gSDij[j] = 0.0f;
			para->getParH(level)->gDxvx[j] = 0.0f;
			para->getParH(level)->gDyvx[j] = 0.0f;
			para->getParH(level)->gDzvx[j] = 0.0f;
			para->getParH(level)->gDxvy[j] = 0.0f;
			para->getParH(level)->gDyvy[j] = 0.0f;
			para->getParH(level)->gDzvy[j] = 0.0f;
			para->getParH(level)->gDxvz[j] = 0.0f;
			para->getParH(level)->gDyvz[j] = 0.0f;
			para->getParH(level)->gDzvz[j] = 0.0f;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	////Taylor Green Vortex uniform
	////like MG
 //   //real uAdvect = real (1. / 250.); //32 nodes -> 250; 40 nodes -> 200; 64 nodes -> 500; 128 nodes -> 1000; 256 nodes -> 2000; 512 nodes -> 4000
 //   const real uAdvect = -0.0016; //32 nodes -> 0.032; 64 nodes -> 0.016; 128 nodes -> 0.008; 256 nodes -> 0.004; 512 nodes -> 0.002

 //   for (int j = 1; j <= numberOfNodes; j++)
 //   {
 //       const real coordX = para->getParH(level)->coordX_SP[j];
 //       const real coordZ = para->getParH(level)->coordZ_SP[j];
 //       const real velocity = para->getVelocity();

 //       para->getParH(level)->rho_SP[j] = real((velocity * velocity) * 3.0 / 4.0 * (cos(coordX * 4.0*PI / gridX) + cos(coordZ * 4.0*PI / gridZ))) * gridZ / gridX;

 //       para->getParH(level)->vy_SP[j] = real(0.0);
 //       para->getParH(level)->vx_SP[j] = real( velocity * sin(coordX * 2.0*PI / gridX) * cos(coordZ * 2.0*PI / gridZ)) + uAdvect * (1.0 + para->getParH(level)->rho_SP[j]);
 //       para->getParH(level)->vz_SP[j] = real(-velocity * cos(coordX * 2.0*PI / gridX) * sin(coordZ * 2.0*PI / gridZ)); // *(real)(gridZ) / (real)(gridX);

 //      //para->getParH(level)->vx_SP[j] = 0.001;
 //      //para->getParH(level)->vy_SP[j] = 0.0;
 //      //para->getParH(level)->vz_SP[j] = 0.001;
 //      //para->getParH(level)->rho_SP[j] = 0.0;
 //      //para->getParH(level)->press_SP[j] = 0.0;

 //       if (para->getCalcMedian()) {
 //           para->getParH(level)->vx_SP_Med[j] = 0.0f;
 //           para->getParH(level)->vy_SP_Med[j] = 0.0f;
 //           para->getParH(level)->vz_SP_Med[j] = 0.0f;
 //           para->getParH(level)->rho_SP_Med[j] = 0.0f;
 //           para->getParH(level)->press_SP_Med[j] = 0.0f;
 //       }
 //       if (para->getUseWale()) {
 //           para->getParH(level)->turbViscosity[j] = 0.0f;
 //           //Debug
 //           para->getParH(level)->gSij[j] = 0.0f;
 //           para->getParH(level)->gSDij[j] = 0.0f;
 //           para->getParH(level)->gDxvx[j] = 0.0f;
 //           para->getParH(level)->gDyvx[j] = 0.0f;
 //           para->getParH(level)->gDzvx[j] = 0.0f;
 //           para->getParH(level)->gDxvy[j] = 0.0f;
 //           para->getParH(level)->gDyvy[j] = 0.0f;
 //           para->getParH(level)->gDzvy[j] = 0.0f;
 //           para->getParH(level)->gDxvz[j] = 0.0f;
 //           para->getParH(level)->gDyvz[j] = 0.0f;
 //           para->getParH(level)->gDzvz[j] = 0.0f;
 //       }


 //       // initial condition

 //       para->getParH(level)->rho_SP[j] = 0.0;

 //       para->getParH(level)->vx_SP[j] = 0.0;
 //       para->getParH(level)->vy_SP[j] = 0.0;
 //       para->getParH(level)->vz_SP[j] = 0.0;

 //   }
	////////////////////////////////////////////////////////////////////////////////
	//2D parabolic test for turbulent channel flow
	//const real uBar = para->getVelocity();				// Bulk velocity computed from DNS results of Kim
	//const real h = 0.5 * gridY;						// half channel width

	//for (int j = 0; j <= numberOfNodes; j++)
	//{
	//	//para->getParH(level)->rho_SP[j] = 
	//	//	(doubflo)((para->getVelocity()*para->getVelocity())*3.0 / 4.0*(cos(para->getParH(level)->coordX_SP[j]*4.0*PI / (doubflo)gridX) + 
	//	//		cos(para->getParH(i)->coordZ_SP[j] *4.0*PI /(real)gridZ)))*(real)(gridZ) / (real)(gridX);
	//	para->getParH(level)->rho_SP[j] =
	//		((real)((para->getVelocity()*para->getVelocity()) * 27.0 *
	//		(cos(para->getParH(level)->coordX_SP[j] * 4.0*PI / (real)gridX) + cos(para->getParH(level)->coordY_SP[j] * 4.0*PI / (real)gridY))) *
	//			(real)(gridY) / (real)(gridX)) +
	//			((real)((para->getVelocity()*para->getVelocity()) * 27.0 *
	//		(cos(para->getParH(level)->coordX_SP[j] * 4.0*PI / (real)gridX) + cos(para->getParH(level)->coordZ_SP[j] * 4.0*PI / (real)gridZ))) *
	//				(real)(gridZ) / (real)(gridX));

	//	//para->getParH(level)->vx_SP[j] = 3.0 * uBar*((para->getParH(level)->coordY_SP[j] / h) - 0.5*((pow(para->getParH(level)->coordY_SP[j], 2.0) / h)));
	//	para->getParH(level)->vx_SP[j] =
	//		3.0 * uBar * (((para->getParH(level)->coordY_SP[j]) / h) - 0.5 / h * ((pow((para->getParH(level)->coordY_SP[j]), 2.0) / h)));

	//	//para->getParH(level)->vy_SP[j] = (real)0.0;
	//	para->getParH(level)->vy_SP[j] =
	//		(real)(-para->getVelocity()*cos((para->getParH(level)->coordX_SP[j] * 2.0*PI / (real)gridX))*sin(para->getParH(level)->coordY_SP[j] * 2.0*PI / (real)gridY));

	//	para->getParH(level)->vz_SP[j] = (real)0.0;
	//	//////////////////////////////////////////////////////////////////////////
	//	if (para->getCalcMedian()) {
	//		para->getParH(level)->vx_SP_Med[j] = 0.0f;
	//		para->getParH(level)->vy_SP_Med[j] = 0.0f;
	//		para->getParH(level)->vz_SP_Med[j] = 0.0f;
	//		para->getParH(level)->rho_SP_Med[j] = 0.0f;
	//		para->getParH(level)->press_SP_Med[j] = 0.0f;
	//	}
	//	if (para->getUseWale()) {
	//		para->getParH(level)->turbViscosity[j] = 0.0f;
	//		//Debug
	//		para->getParH(level)->gSij[j] = 0.0f;
	//		para->getParH(level)->gSDij[j] = 0.0f;
	//		para->getParH(level)->gDxvx[j] = 0.0f;
	//		para->getParH(level)->gDyvx[j] = 0.0f;
	//		para->getParH(level)->gDzvx[j] = 0.0f;
	//		para->getParH(level)->gDxvy[j] = 0.0f;
	//		para->getParH(level)->gDyvy[j] = 0.0f;
	//		para->getParH(level)->gDzvy[j] = 0.0f;
	//		para->getParH(level)->gDxvz[j] = 0.0f;
	//		para->getParH(level)->gDyvz[j] = 0.0f;
	//		para->getParH(level)->gDzvz[j] = 0.0f;
	//	}
	//}
	////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	// Lenz: 3D- Taylor-Green Vortex
	//////////////////////////////////////////////////////////////////////////
	//for (int j = 1; j <= numberOfNodes; j++)
	//{
	//    const real coordX = para->getParH(level)->coordX_SP[j];
	//    const real coordY = para->getParH(level)->coordY_SP[j];
	//    const real coordZ = para->getParH(level)->coordZ_SP[j];
	//    const real velocity = para->getVelocity();

	//    real A =  1.0;
	//    real B =  1.0;
	//    real C = -2.0;
	//    real a = 2.0 * PI;
	//    real b = 2.0 * PI;
	//    real c = 2.0 * PI;

	//    //para->getParH(level)->rho_SP[j] = 0.0;

	//    //para->getParH(level)->vx_SP[j] = 0.05 * A * cos( a * coordX ) * sin( b * coordY ) * sin( c * coordZ ) + 0.05;
	//    //para->getParH(level)->vy_SP[j] = 0.05 * B * sin( a * coordX ) * cos( b * coordY ) * sin( c * coordZ );
	//    //para->getParH(level)->vz_SP[j] = 0.05 * C * sin( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );

	//    para->getParH(level)->rho_SP[j] = 0.0;//1.0 / 16.0 * ( cos( 2.0 * a * coordX ) + cos( 2.0 * b * coordY ) ) * ( cos( 2.0 * c * coordZ ) + 2.0 );

	//    para->getParH(level)->vx_SP[j] =   0.05 * sin( a * coordX ) * cos( b * coordY ) * cos( c * coordZ );
	//    para->getParH(level)->vy_SP[j] = - 0.05 * cos( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );
	//    para->getParH(level)->vz_SP[j] =   0.0;

	//   //para->getParH(level)->vx_SP[j] = 0.001;
	//   //para->getParH(level)->vy_SP[j] = 0.0;
	//   //para->getParH(level)->vz_SP[j] = 0.001;
	//   //para->getParH(level)->rho_SP[j] = 0.0;
	//   //para->getParH(level)->press_SP[j] = 0.0;

	//    if (para->getCalcMedian()) {
	//        para->getParH(level)->vx_SP_Med[j] = 0.0f;
	//        para->getParH(level)->vy_SP_Med[j] = 0.0f;
	//        para->getParH(level)->vz_SP_Med[j] = 0.0f;
	//        para->getParH(level)->rho_SP_Med[j] = 0.0f;
	//        para->getParH(level)->press_SP_Med[j] = 0.0f;
	//    }
	//    if (para->getUseWale()) {
	//        para->getParH(level)->turbViscosity[j] = 0.0f;
	//        //Debug
	//        para->getParH(level)->gSij[j] = 0.0f;
	//        para->getParH(level)->gSDij[j] = 0.0f;
	//        para->getParH(level)->gDxvx[j] = 0.0f;
	//        para->getParH(level)->gDyvx[j] = 0.0f;
	//        para->getParH(level)->gDzvx[j] = 0.0f;
	//        para->getParH(level)->gDxvy[j] = 0.0f;
	//        para->getParH(level)->gDyvy[j] = 0.0f;
	//        para->getParH(level)->gDzvy[j] = 0.0f;
	//        para->getParH(level)->gDxvz[j] = 0.0f;
	//        para->getParH(level)->gDyvz[j] = 0.0f;
	//        para->getParH(level)->gDzvz[j] = 0.0f;
	//    }

	//}
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	// 3D-Taylor-Green Vortex from Paper:"Dynamical selective filtering for the Lattice Boltzmann Method" (Simon Marié, Xavier Gloerfelt)
	//////////////////////////////////////////////////////////////////////////
	//for (int j = 1; j <= numberOfNodes; j++)
	//{
	//    const real coordX = para->getParH(level)->coordX_SP[j];
	//    const real coordY = para->getParH(level)->coordY_SP[j];
	//    const real coordZ = para->getParH(level)->coordZ_SP[j];
	//    const real velocity = para->getVelocity();

	//    real A = 2.0;//1.0;//
	//	real B = 2.0;//1.0;//
	//	real C = 2.0;//1.0;//
	//	real a = 1.0;//A * PI;//1.0;//
	//	real b = 1.0;//B * PI;//1.0;//
	//	real c = 1.0;//C * PI;//1.0;//
	//	//////////////////////////////////////////////////////////////////////////

	//    para->getParH(level)->rho_SP[j] = /*1.0 +*/ 3.0 * ((velocity * velocity) / 16.0 * ( cos( 2.0 * a * coordX ) + cos( 2.0 * b * coordY ) ) * ( cos( 2.0 * c * coordZ ) + 2.0 ) );

	//    para->getParH(level)->vx_SP[j] =  velocity * sin( a * coordX ) * cos( b * coordY ) * cos( c * coordZ );
	//    para->getParH(level)->vy_SP[j] = -velocity * cos( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );
	//    para->getParH(level)->vz_SP[j] =  0.0; //0.00001 * velocity * sin(4.0 * a * coordX) * cos(4.0 * b * coordY) * cos(4.0 * c * coordZ); //0.0;

	//    if (para->getCalcMedian()) {
	//        para->getParH(level)->vx_SP_Med[j] = 0.0f;
	//        para->getParH(level)->vy_SP_Med[j] = 0.0f;
	//        para->getParH(level)->vz_SP_Med[j] = 0.0f;
	//        para->getParH(level)->rho_SP_Med[j] = 0.0f;
	//        para->getParH(level)->press_SP_Med[j] = 0.0f;
	//    }
	//    if (para->getUseWale()) {
	//        para->getParH(level)->turbViscosity[j] = 0.0f;
	//        //Debug
	//        para->getParH(level)->gSij[j] = 0.0f;
	//        para->getParH(level)->gSDij[j] = 0.0f;
	//        para->getParH(level)->gDxvx[j] = 0.0f;
	//        para->getParH(level)->gDyvx[j] = 0.0f;
	//        para->getParH(level)->gDzvx[j] = 0.0f;
	//        para->getParH(level)->gDxvy[j] = 0.0f;
	//        para->getParH(level)->gDyvy[j] = 0.0f;
	//        para->getParH(level)->gDzvy[j] = 0.0f;
	//        para->getParH(level)->gDxvz[j] = 0.0f;
	//        para->getParH(level)->gDyvz[j] = 0.0f;
	//        para->getParH(level)->gDzvz[j] = 0.0f;
	//    }

	//}
	//////////////////////////////////////////////////////////////////////////
}


void GridProvider::setPressSizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->QPress.kQ = sizePerLevel;
	para->getParD(level)->QPress.kQ = sizePerLevel;
	para->getParH(level)->kPressQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kPressQread = sizePerLevel * para->getD3Qxx();
}


void GridProvider::setVelocitySizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->Qinflow.kQ = sizePerLevel;
	para->getParD(level)->Qinflow.kQ = sizePerLevel;
	para->getParH(level)->kInflowQ = sizePerLevel;
	para->getParD(level)->kInflowQ = sizePerLevel;
	para->getParH(level)->kInflowQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kInflowQread = sizePerLevel * para->getD3Qxx();
}

void GridProvider::setOutflowSizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->Qoutflow.kQ = sizePerLevel;
	para->getParD(level)->Qoutflow.kQ = sizePerLevel;
	para->getParH(level)->kOutflowQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kOutflowQread = sizePerLevel * para->getD3Qxx();
}

void GridProvider::allocAndCopyForcing()
{
	cudaMemoryManager->cudaAllocForcing();
	cudaMemoryManager->cudaCopyForcingToDevice();
}

void GridProvider::freeMemoryOnHost()
{
	for (int level = para->getCoarse(); level <= para->getFine(); level++)
	{
		cudaMemoryManager->cudaFreeCoord(level);
		cudaMemoryManager->cudaFreeSP(level);
	}
}

void GridProvider::cudaCopyDataToHost(int level)
{
	cudaMemoryManager->cudaCopyDataToHost(level);
}
