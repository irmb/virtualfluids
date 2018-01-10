#include "Init/DefineGrid.h"
#include "Init/ReadGeometry.h"
#include "Temperature/FindTemperature.h"
#include "FindInterface/FindInterface.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

////////////////////////////////////////////////////////////////////////////////
void defineGrid(Parameter* para, Communicator* comm)
{
	for (int lev=para->getFine(); lev >= para->getCoarse(); lev--)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////
		// Allocate Host Memory
		///////////////////////////////////////////////////////////////////////////////////////////////////
		para->cudaAllocFull(lev);
		///////////////////////////////////////////////////////////////////////////////////////////////////
		if (para->getDiffOn()==true)
		{
			checkCudaErrors( cudaMallocHost((void**) &(para->getParH(lev)->Conc_Full ), para->getParH(lev)->mem_size_real));
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////
		if (lev==para->getCoarse())
		{
			///////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getReadGeo()==true)
			{
				std::cout << "read geometry...\n" ;
				readGeometry(para, comm, lev, para->getGeometryFileC());
				std::cout << "done.\n";
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////
			for (unsigned int k=0; k<para->getParH(lev)->nz; k++)
			{
				for (unsigned int j=0; j<para->getParH(lev)->ny; j++)
				{
					for (unsigned int i=0; i<para->getParH(lev)->nx; i++)
					{
						int m = para->getParH(lev)->nx*(para->getParH(lev)->ny*k + j) + i;
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						if (para->getReadGeo()==true)
						{
							if (  i <= STARTOFFX || i >= para->getParH(lev)->gridNX + STARTOFFX - 1 
								|| j <= STARTOFFY || j >= para->getParH(lev)->gridNY + STARTOFFY - 1 
								|| k <= STARTOFFZ || k >= para->getParH(lev)->gridNZ + STARTOFFZ - 1 )
							{
								para->getParH(lev)->geo[m] = GEO_VOID;
							}
							//else if (  (k >= STARTOFFZ +1) && (k <= STARTOFFZ + 20)  )
							//{
							//   para->getParH(lev)->geo[m] = GEO_FLUID;
							//}
						}
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						else
						{
							////////////////////////////////////////////////////////////////////////////////////////////////
							//Test
							////unsigned int centerX = para->getParH(lev)->gridNX / 2 + STARTOFFX;
							////unsigned int centerY = para->getParH(lev)->gridNY / 2 + STARTOFFY;
							////unsigned int centerZ = para->getParH(lev)->gridNZ / 4 + STARTOFFZ;
							//real centerX = para->getParH(lev)->gridNX / 2. + STARTOFFX - 0.5;
							//real centerY = para->getParH(lev)->gridNY / 2. + STARTOFFY - 0.5;
							//real centerZ = para->getParH(lev)->gridNZ / 4. + STARTOFFZ - 0.5;
							//real      radius  = para->getParH(lev)->gridNY / 10.;//2.56f;
							//////unsigned int distSq  = (centerX-i)*(centerX-i)+(centerY-j)*(centerY-j);
							////unsigned int distSq  = (centerX-i)*(centerX-i)+(centerY-j)*(centerY-j)+(centerZ-k)*(centerZ-k);
							//real distSq  = (centerX-i)*(centerX-i)+(centerY-j)*(centerY-j)+(centerZ-k)*(centerZ-k);
							//real      radiSq  = radius*radius;
							////////////////////////////////////////////////////////////////////////////////////////////////

							if (   i < STARTOFFX || i > para->getParH(lev)->gridNX + STARTOFFX - 1 
								|| j < STARTOFFY || j > para->getParH(lev)->gridNY + STARTOFFY - 1 
								|| k < STARTOFFZ || k > para->getParH(lev)->gridNZ + STARTOFFZ - 1 )
							{
								para->getParH(lev)->geo[m] = GEO_VOID;
								if (para->getDiffOn()==true)
								{
									para->getParH(lev)->Conc_Full[m] = 0.0;
								}
							}
							//else if (  j == STARTOFFY  || j == para->getParH(lev)->gridNY + STARTOFFY - 1 )
							//{
							//   para->getParH(lev)->geo[m] = GEO_SOLID;
							//   if (ic.diffOn==true)
							//   {
							//      para->getParH(lev)->Conc_Full[m] = 0.0;
							//   }
							//}
							//else if (  k == STARTOFFZ || k == para->getParH(lev)->gridNZ + STARTOFFZ - 1)
							//{
							//   para->getParH(lev)->geo[m] = GEO_SOLID;
							//   if (ic.diffOn==true)
							//   {
							//      para->getParH(lev)->Conc_Full[m] = 0.0;
							//   }
							//}
							//else if (  i == STARTOFFX+1  || i == para->getParH(lev)->gridNX + STARTOFFX - 2 )
							//{
							//   para->getParH(lev)->geo[m] = GEO_VELO;
							//   if (ic.diffOn==true)
							//   {
							//      para->getParH(lev)->Conc_Full[m] = 0.0;
							//   }
							//}
							//else if ((i >= para->getParH(lev+1)->XdistKn + STARTOFFX + 3) && (i <= (para->getParH(lev+1)->XdistKn + para->getParH(lev+1)->gridNX / 2) + STARTOFFX - 3) &&
							//         (j >= para->getParH(lev+1)->YdistKn + STARTOFFY + 3) && (j <= (para->getParH(lev+1)->YdistKn + para->getParH(lev+1)->gridNY / 2) + STARTOFFY - 3) &&
							//         (k >= para->getParH(lev+1)->ZdistKn + STARTOFFZ /*+ 3*/) && (k <= (para->getParH(lev+1)->ZdistKn + para->getParH(lev+1)->gridNZ / 2) + STARTOFFZ - 3) )
							//{
							//   para->getParH(lev)->geo[m] = GEO_VOID;
							//   if (para->getDiffOn()==true)
							//   {
							//      para->getParH(lev)->Conc_Full[m] = 0.0;
							//   }
							//}
							//else if ((i >= para->getParH(lev)->gridNX / 4 + STARTOFFX + 4) && (i <= (para->getParH(lev)->gridNX * 3 / 4) + STARTOFFX - 4) &&
							//         (j >= para->getParH(lev)->gridNY / 4 + STARTOFFY + 4) && (j <= (para->getParH(lev)->gridNY * 3 / 4) + STARTOFFY - 4) &&
							//         (k >= para->getParH(lev)->gridNZ / 8 + STARTOFFZ + 4) && (k <= (para->getParH(lev)->gridNZ * 5 / 8) + STARTOFFZ - 4) && maxlevel>1)
							//{
							//   para->getParH(lev)->geo[m] = GEO_SOLID;
							//   if (ic.diffOn==true)
							//   {
							//      para->getParH(lev)->Conc_Full[m] = 0.0;
							//   }
							//}
							//else if (distSq < /*>*/ radiSq)
							//{
							//   para->getParH(lev)->geo[m] = GEO_SOLID;
							//      if (ic.diffOn==true)
							//      {
							//         para->getParH(lev)->Conc_Full[m] = 0.0;
							//      }
							//}
							//else if (  i <= STARTOFFX + 30 && i >=STARTOFFX + 20 
							//        && j <= STARTOFFY + 30 && j >=STARTOFFY + 20 
							//        && k <= STARTOFFZ + 30 && k >=STARTOFFZ + 20 )
							//{
							//   para->getParH(lev)->geo[m] = GEO_FLUID;
							//   para->getParH(lev)->Conc_Full[m] = 1.0;
							//}
							//else if (   i <= STARTOFFX + 30 && i >=STARTOFFX + 20 
							//         && j <= STARTOFFY + 30 && j >=STARTOFFY + 20 
							//         && k <= STARTOFFZ + 30 && k >=STARTOFFZ + 20 )
							//{
							//   para->getParH(lev)->geo[m] = GEO_FLUID;
							//   para->getParH(lev)->Conc_Full[m] = 1.0;
							//}
							else
							{
								para->getParH(lev)->geo[m] = GEO_FLUID;
								if (para->getDiffOn()==true)
								{
									para->getParH(lev)->Conc_Full[m] = (real)para->getTemperatureInit();
								}
							}
							//if (i == STARTOFFZ)
							//{
							// para->getParH(lev)->geo[m] = GEO_SOLID;
							// if (para->getDiffOn()==true)
							// {
							//	 para->getParH(lev)->Conc_Full[m] = 0.0;
							// }
							//}

						}
						para->getParH(lev)->k[m] = 0;
					}
				}
			}
		}
		else if (lev==para->getFine() && para->getMaxLevel()>=1)
		{
			for (unsigned int k=0; k<para->getParH(lev)->nz; k++)
			{
				for (unsigned int j=0; j<para->getParH(lev)->ny; j++)
				{
					for (unsigned int i=0; i<para->getParH(lev)->nx; i++)
					{
						//unsigned int centerX = para->getParH(lev)->gridNX / 2 + STARTOFFX;
						//unsigned int centerY = para->getParH(lev)->gridNY / 2 + STARTOFFY;
						//unsigned int centerZ = para->getParH(lev)->gridNZ / 2 + STARTOFFZ;
						//real        radius  = para->getParH(lev)->gridNY / 5.f;//2.56f;
						//real centerX = para->getParH(lev)->gridNX / 2.f + STARTOFFX - 0.5f;
						//real centerY = para->getParH(lev)->gridNY / 2.f + STARTOFFY - 0.5f;
						//real centerZ = para->getParH(lev)->gridNZ / 2.f + STARTOFFZ - 0.5f;
						//real      radius  = para->getParH(lev)->gridNY / 5.f;//2.56f;

						int m = para->getParH(lev)->nx*(para->getParH(lev)->ny*k + j) + i;

						////unsigned int distSq = (centerX-i)*(centerX-i)+(centerY-j)*(centerY-j)+(centerZ-k)*(centerZ-k);
						//real distSq = (centerX-i)*(centerX-i)+(centerY-j)*(centerY-j)+(centerZ-k)*(centerZ-k);
						//real radiSq = radius*radius;

						//diff stuff
						//real mradsq = (real)((real)i-(STARTOFFX + 30)) * (real)((real)i-(STARTOFFX + 30)) + (real)((real)j-(STARTOFFY + 30)) *  (real)((real)j-(STARTOFFY + 30)) +(real) ((real)k-(STARTOFFZ + 30)) * (real) ((real)k-(STARTOFFZ + 30)); 

						if (  i <  STARTOFFX || i >  para->getParH(lev)->gridNX + STARTOFFX - 1 
							|| j <  STARTOFFY || j >  para->getParH(lev)->gridNY + STARTOFFY - 1 
							|| k <  STARTOFFZ || k >  para->getParH(lev)->gridNZ + STARTOFFZ - 1 )
						{
							para->getParH(lev)->geo[m] = GEO_VOID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						//else if (i = STARTOFFX)
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						//else if ((i >  STARTOFFX )     && (i <= para->getParH(lev)->gridNX + STARTOFFX - 38) &&
						//         (j >= 19 + STARTOFFY) && (j <= para->getParH(lev)->gridNY + STARTOFFY - 19) &&
						//         (k >= 30 + STARTOFFZ) && (k <= para->getParH(lev)->gridNZ + STARTOFFZ - 30) )
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						//else if ((i >= STARTOFFX )                                 && (i <= (para->getParH(lev)->gridNX * 4 / 8) + STARTOFFX) &&
						//         (j >= para->getParH(lev)->gridNY     /  4 + STARTOFFY) && (j <= (para->getParH(lev)->gridNY *  3 /  4) + STARTOFFY) &&
						//         (k >= para->getParH(lev)->gridNZ * 8 / 20 + STARTOFFZ) && (k <= (para->getParH(lev)->gridNZ * 12 / 20) + STARTOFFZ) )
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						////diff stuff
						//else if (  mradsq < 100.f )
						//{
						//   para->getParH(lev)->geo[m] = GEO_FLUID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 1.0f - mradsq * 0.01f;
						//   }
						//}
						//else if (  i <= STARTOFFX + 30 && i >=STARTOFFX + 20 
						//        && j <= STARTOFFY + 30 && j >=STARTOFFY + 20 
						//        && k <= STARTOFFZ + 30 && k >=STARTOFFZ + 20 )
						//{
						//   para->getParH(lev)->geo[m] = GEO_FLUID;
						//   para->getParH(lev)->Conc_Full[m] = 1.0;
						//}
						//else if (distSq < radiSq)
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						else
						{
							para->getParH(lev)->geo[m] = GEO_FLUID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						//if (i == STARTOFFX)
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						if (k == STARTOFFZ)
						{
							para->getParH(lev)->geo[m] = GEO_SOLID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}


						para->getParH(lev)->k[m] = 0;
					}
				}
			}
		}
		else if( lev > para->getCoarse() && lev < para->getFine() )
		{
			for (unsigned int k=0; k<para->getParH(lev)->nz; k++)
			{
				for (unsigned int j=0; j<para->getParH(lev)->ny; j++)
				{
					for (unsigned int i=0; i<para->getParH(lev)->nx; i++)
					{
						int m = para->getParH(lev)->nx*(para->getParH(lev)->ny*k + j) + i;

						if (  i < STARTOFFX || i > para->getParH(lev)->gridNX + STARTOFFX - 1 
							|| j < STARTOFFY || j > para->getParH(lev)->gridNY + STARTOFFY - 1 
							|| k < STARTOFFZ || k > para->getParH(lev)->gridNZ + STARTOFFZ - 1 )
						{
							para->getParH(lev)->geo[m] = GEO_VOID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						else if ((i >= para->getParH(lev+1)->XdistKn + STARTOFFX + 3) && (i <= (para->getParH(lev+1)->XdistKn + para->getParH(lev+1)->gridNX / 2) + STARTOFFX - 3) &&
							(j >= para->getParH(lev+1)->YdistKn + STARTOFFY + 3) && (j <= (para->getParH(lev+1)->YdistKn + para->getParH(lev+1)->gridNY / 2) + STARTOFFY - 3) &&
							(k >= para->getParH(lev+1)->ZdistKn + STARTOFFZ + 3) && (k <= (para->getParH(lev+1)->ZdistKn + para->getParH(lev+1)->gridNZ / 2) + STARTOFFZ - 3) )
						{
							para->getParH(lev)->geo[m] = GEO_VOID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						//else if (i = STARTOFFX)
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//   if (ic.diffOn==true)
						//   {
						//      para->getParH(lev)->Conc_Full[m] = 0.0;
						//   }
						//}
						//else if ((i >= para->getParH(lev)->gridNX / 4 + STARTOFFX + 3) && (i <= (para->getParH(lev)->gridNX * 3 / 4) + STARTOFFX - 3) &&
						//         (j >= para->getParH(lev)->gridNY / 4 + STARTOFFY + 3) && (j <= (para->getParH(lev)->gridNY * 3 / 4) + STARTOFFY - 3) &&
						//         (k >= para->getParH(lev)->gridNZ / 8 + STARTOFFZ + 3) && (k <= (para->getParH(lev)->gridNZ * 5 / 8) + STARTOFFZ - 3) )
						//{
						//   para->getParH(lev)->geo[m] = GEO_SOLID;
						//}
						else
						{
							para->getParH(lev)->geo[m] = GEO_FLUID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						if (i == STARTOFFX)
						{
							para->getParH(lev)->geo[m] = GEO_SOLID;
							if (para->getDiffOn()==true)
							{
								para->getParH(lev)->Conc_Full[m] = 0.0;
							}
						}
						para->getParH(lev)->k[m] = 0;
					}
				}
			}
		}

		//std::cout << "read geoFull..." ;
		//readVFgeoFull(para, "D:/temp/gpuBenchmarkCylinder/GPU/geoFull.dat");
		//std::cout << "done.\n";

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Buffer GEO
		//geo_sbuf_t.setSize(1,para->getParH(0)->sizePlaneXY);
		//geo_rbuf_t.setSize(1,para->getParH(0)->sizePlaneXY);
		//geo_sbuf_b.setSize(1,para->getParH(0)->sizePlaneXY);
		//geo_rbuf_b.setSize(1,para->getParH(0)->sizePlaneXY);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Exchange GEO
		//if (numprocs>1)
		//{
		//   exchangeDataGeo(lev);
		//}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->setSizeMatSparse(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->cudaAllocSP(lev);
		//F3
		para->cudaAllocF3SP(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Buffer
		//sbuf_t.setSize(27,para->getParH(0)->sizePlaneST);
		//rbuf_t.setSize(27,para->getParH(0)->sizePlaneRT);
		//sbuf_b.setSize(27,para->getParH(0)->sizePlaneSB);
		//rbuf_b.setSize(27,para->getParH(0)->sizePlaneRB);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (para->getDiffOn()==true)
		{
			para->cudaAllocConc(lev);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->fillSparse(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->cudaCopySP(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (para->getDiffOn()==true)
		{
			std::cout << "Maikes Thermo-Stuff...\n" ;
			initTemperatur(para, lev);//thermostuff(lev);
			std::cout << "done.\n";
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if ( lev < para->getFine() && para->getMaxLevel()>=1)
		{
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocInterfaceCF(lev);
			para->cudaAllocInterfaceFC(lev);
			para->cudaAllocInterfaceOffCF(lev);
			para->cudaAllocInterfaceOffFC(lev);
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//Find Interpolation Interface
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout << "Anzahl der CF-Interpolationzellen vorher:" <<  para->getParH(lev)->K_CF << "\n";
			std::cout << "Anzahl der FC-Interpolationzellen vorher:" <<  para->getParH(lev)->K_FC << "\n";
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			interpolation(para->getParH(lev)->intCF,     para->getParH(lev)->intFC, 
				para->getParH(lev)->gridNX,    para->getParH(lev)->gridNY,    para->getParH(lev)->gridNZ, 
				para->getParH(lev+1)->gridNX,  para->getParH(lev+1)->gridNY,  para->getParH(lev+1)->gridNZ, 
				para->getParH(lev+1)->XdistKn, para->getParH(lev+1)->YdistKn, para->getParH(lev+1)->ZdistKn,
				para->getParH(lev)->k,         para->getParH(lev+1)->k,       para->getParH(lev)->need_interface,
				para->getParH(lev)->offCF,     para->getParH(lev)->offFC);
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout << "Anzahl der CF-Interpolationzellen nachher:" <<  para->getParH(lev)->intCF.kCF << "\n";
			std::cout << "Anzahl der FC-Interpolationzellen nachher:" <<  para->getParH(lev)->intFC.kFC << "\n";
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//Test
			//for (int test = 0; test < para->getParH(lev)->intCF.kCF; test++)
			//{
			// para->getParH(lev)->offCF.xOffCF[test] = -para->getParH(lev)->offCF.xOffCF[test];
			// para->getParH(lev)->offCF.yOffCF[test] = -para->getParH(lev)->offCF.yOffCF[test];
			// para->getParH(lev)->offCF.zOffCF[test] = -para->getParH(lev)->offCF.zOffCF[test];
			//}
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyInterfaceCF(lev);
			para->cudaCopyInterfaceFC(lev);
			para->cudaCopyInterfaceOffCF(lev);
			para->cudaCopyInterfaceOffFC(lev);
			//////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//if (para->getMyID() == para->getPressInID())       setSizeOfPlane(para, 0, para->getPressInZ());
	//else if(para->getMyID() == para->getPressOutID())  setSizeOfPlane(para, 0, para->getPressOutZ());
	////////////////////////////////////////////////////////////////////////////////////////////////////////

}
