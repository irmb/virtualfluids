#include "Parameter.h"

//#include <cuda_runtime.h>
//#include <helper_cuda.h>
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#ifdef WIN32
//   #include <Winsock2.h>
//#endif
//lib for windows Ws2_32.lib

SPtr<Parameter> Parameter::make()
{
    return SPtr<Parameter>(new Parameter());
}

Parameter::Parameter()
{
}
Parameter* Parameter::instanz = 0;
Parameter* Parameter::getInstanz()
{
	if( instanz == 0 )
		instanz = new Parameter();
	return instanz;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//init-method
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::initParameter()
{
	factor_gridNZ  = 2;
	coarse         = 0;
	fine           = this->maxlevel;
	parH.resize(this->maxlevel+1);
	parD.resize(this->maxlevel+1);

	//host
	for (int i = coarse; i <= fine; i++)
	{
		parH[i]                        = new ParameterStruct;
		parH[i]->numberofthreads       = 64;// 128;
		parH[i]->gridNX                = getGridX().at(i);
		parH[i]->gridNY                = getGridY().at(i);
		parH[i]->gridNZ                = getGridZ().at(i);
		parH[i]->vis                   = ic.vis*pow(2.f,i);
		parH[i]->diffusivity           = ic.Diffusivity*pow(2.f,i);
		parH[i]->omega                 = 1.0f/(3.0f*parH[i]->vis+0.5f);//omega :-) not s9 = -1.0f/(3.0f*parH[i]->vis+0.5f);//
		parH[i]->nx                    = parH[i]->gridNX + 2 * STARTOFFX;
		parH[i]->ny                    = parH[i]->gridNY + 2 * STARTOFFY;
		parH[i]->nz                    = parH[i]->gridNZ + 2 * STARTOFFZ;
		parH[i]->size_Mat              = parH[i]->nx * parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXY           = parH[i]->nx * parH[i]->ny;
		parH[i]->sizePlaneYZ           = parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXZ           = parH[i]->nx * parH[i]->nz;
		parH[i]->mem_size_real      = sizeof(real     ) * parH[i]->size_Mat;
		parH[i]->mem_size_int          = sizeof(unsigned int) * parH[i]->size_Mat;
		parH[i]->mem_size_bool         = sizeof(bool        ) * parH[i]->size_Mat;
		parH[i]->mem_size_real_yz   = sizeof(real     ) * parH[i]->ny * parH[i]->nz;
		parH[i]->evenOrOdd             = true;
		parH[i]->startz                = parH[i]->gridNZ * ic.myid;
		parH[i]->endz                  = parH[i]->gridNZ * ic.myid + parH[i]->gridNZ;
		parH[i]->Lx                    = (real)((1.f*parH[i]->gridNX - 1.f)/(pow(2.f,i)));
		parH[i]->Ly                    = (real)((1.f*parH[i]->gridNY - 1.f)/(pow(2.f,i)));
		parH[i]->Lz                    = (real)((1.f*parH[i]->gridNZ - 1.f)/(pow(2.f,i)));
		parH[i]->dx                    = (real)(1.f/(pow(2.f,i)));
		parH[i]->XdistKn               = getDistX().at(i);
		parH[i]->YdistKn               = getDistY().at(i);
		parH[i]->ZdistKn               = getDistZ().at(i);
		if (i==coarse)
		{
			parH[i]->distX                 = (real)getDistX().at(i);
			parH[i]->distY                 = (real)getDistY().at(i);
			parH[i]->distZ                 = (real)getDistZ().at(i);
			parH[i]->mTtoWx                = (real)1.0f;
			parH[i]->mTtoWy                = (real)1.0f;
			parH[i]->mTtoWz                = (real)1.0f;
			parH[i]->cTtoWx                = (real)0.0f;
			parH[i]->cTtoWy                = (real)0.0f;
			parH[i]->cTtoWz                = (real)0.0f;
			////MGs Trafo///////////////////////////////////////////////////////////////
			//parH[i]->cStartx               = (real)parH[i]->XdistKn;
			//parH[i]->cStarty               = (real)parH[i]->XdistKn;
			//parH[i]->cStartz               = (real)parH[i]->XdistKn;
			////////////////////////////////////////////////////////////////////////////
		} 
		else
		{
			//Geller
			parH[i]->distX                 = ((real)getDistX().at(i) + 0.25f) * parH[i-1]->dx;
			parH[i]->distY                 = ((real)getDistY().at(i) + 0.25f) * parH[i-1]->dx;
			parH[i]->distZ                 = ((real)getDistZ().at(i) + 0.25f) * parH[i-1]->dx;
			//parH[i]->distX                 = ((real)getDistX().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distX;
			//parH[i]->distY                 = ((real)getDistY().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distY;
			//parH[i]->distZ                 = ((real)getDistZ().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distZ;
			parH[i]->mTtoWx                = (real)pow(0.5f,i);
			parH[i]->mTtoWy                = (real)pow(0.5f,i);
			parH[i]->mTtoWz                = (real)pow(0.5f,i);
			parH[i]->cTtoWx                = (real)(STARTOFFX/2.f + (parH[i]->gridNX+1.f)/4.f); //funzt nur für zwei level
			parH[i]->cTtoWy                = (real)(STARTOFFY/2.f + (parH[i]->gridNY+1.f)/4.f); //funzt nur für zwei level
			parH[i]->cTtoWz                = (real)(STARTOFFZ/2.f + (parH[i]->gridNZ+1.f)/4.f); //funzt nur für zwei level
			////MGs Trafo///////////////////////////////////////////////////////////////
			//parH[i]->cStartx               = (real)parH[i]->XdistKn;
			//parH[i]->cStarty               = (real)parH[i]->XdistKn;
			//parH[i]->cStartz               = (real)parH[i]->XdistKn;
			////////////////////////////////////////////////////////////////////////////
		}
		parH[i]->need_interface[INTERFACE_E]=getNeedInterface().at(INTERFACE_E);
		parH[i]->need_interface[INTERFACE_W]=getNeedInterface().at(INTERFACE_W);
		parH[i]->need_interface[INTERFACE_N]=getNeedInterface().at(INTERFACE_N);
		parH[i]->need_interface[INTERFACE_S]=getNeedInterface().at(INTERFACE_S);
		parH[i]->need_interface[INTERFACE_T]=getNeedInterface().at(INTERFACE_T);
		parH[i]->need_interface[INTERFACE_B]=getNeedInterface().at(INTERFACE_B);
	}

	//device
	for (int i = coarse; i <= fine; i++)
	{
		parD[i]                        = new ParameterStruct;
		parD[i]->numberofthreads       = parH[i]->numberofthreads;
		parD[i]->gridNX                = parH[i]->gridNX;
		parD[i]->gridNY                = parH[i]->gridNY;
		parD[i]->gridNZ                = parH[i]->gridNZ;
		parD[i]->vis                   = parH[i]->vis;
		parD[i]->diffusivity           = parH[i]->diffusivity;
		parD[i]->omega                 = parH[i]->omega;
		parD[i]->nx                    = parH[i]->nx;
		parD[i]->ny                    = parH[i]->ny;
		parD[i]->nz                    = parH[i]->nz;
		parD[i]->size_Mat              = parH[i]->size_Mat;
		parD[i]->sizePlaneXY           = parH[i]->sizePlaneXY;
		parD[i]->sizePlaneYZ           = parH[i]->sizePlaneYZ;
		parD[i]->sizePlaneXZ           = parH[i]->sizePlaneXZ;
		parD[i]->mem_size_real      = sizeof(real     ) * parD[i]->size_Mat;
		parD[i]->mem_size_int          = sizeof(unsigned int) * parD[i]->size_Mat;
		parD[i]->mem_size_bool         = sizeof(bool        ) * parD[i]->size_Mat;
		parD[i]->mem_size_real_yz   = sizeof(real     ) * parD[i]->ny * parD[i]->nz;
		parD[i]->evenOrOdd             = parH[i]->evenOrOdd;
		parD[i]->startz                = parH[i]->startz;
		parD[i]->endz                  = parH[i]->endz;
		parD[i]->Lx                    = parH[i]->Lx;
		parD[i]->Ly                    = parH[i]->Ly;
		parD[i]->Lz                    = parH[i]->Lz;
		parD[i]->dx                    = parH[i]->dx;
		parD[i]->XdistKn               = parH[i]->XdistKn;
		parD[i]->YdistKn               = parH[i]->YdistKn;
		parD[i]->ZdistKn               = parH[i]->ZdistKn;
		parD[i]->distX                 = parH[i]->distX;
		parD[i]->distY                 = parH[i]->distY;
		parD[i]->distZ                 = parH[i]->distZ;
	}

	//Interface
	//comment out for geller
	//for (int i = coarse; i < fine; i++)
	//{
	//   initInterfaceParameter(i);
	//}
}
void Parameter::setSizeMatSparse(int level)
{
	parH[level]->size_Mat_SP = 1;
	parD[level]->size_Mat_SP = 1;
	parH[level]->sizePlaneSB = 0;
	parH[level]->sizePlaneST = 0;
	parH[level]->sizePlaneRB = 0;
	parH[level]->sizePlaneRT = 0;
	parH[level]->isSetSendB  = false;
	parH[level]->isSetSendT  = false;
	parH[level]->isSetRecvB  = false;
	parH[level]->isSetRecvT  = false;
	unsigned int mm[8];

	for (unsigned int k=1; k<parH[level]->gridNZ + 2 * STARTOFFZ - 1; k++)
	{
		for (unsigned int j=1; j<parH[level]->gridNY + 2 * STARTOFFY - 1; j++)
		{
			for (unsigned int i=1; i<parH[level]->gridNX + 2 * STARTOFFX - 1; i++)
			{
				mm[0]= parH[level]->nx*(parH[level]->ny*k + j) + i;
				mm[1]= mm[0]                                                  -1; //W
				mm[2]= mm[0]                                  -parH[level]->nx-1; //SW
				mm[3]= mm[0]                                  -parH[level]->nx;   //S
				mm[4]= mm[0]-(parH[level]->nx*parH[level]->ny);                   //B
				mm[5]= mm[0]-(parH[level]->nx*parH[level]->ny)                -1; //BW
				mm[6]= mm[0]-(parH[level]->nx*parH[level]->ny)-parH[level]->nx;   //BS
				mm[7]= mm[0]-(parH[level]->nx*parH[level]->ny)-parH[level]->nx-1; //BSW

				if ( parH[level]->geo[mm[0]] != GEO_VOID ||
					parH[level]->geo[mm[1]] != GEO_VOID ||
					parH[level]->geo[mm[2]] != GEO_VOID ||
					parH[level]->geo[mm[3]] != GEO_VOID ||
					parH[level]->geo[mm[4]] != GEO_VOID ||
					parH[level]->geo[mm[5]] != GEO_VOID ||
					parH[level]->geo[mm[6]] != GEO_VOID ||
					parH[level]->geo[mm[7]] != GEO_VOID )
				{
					//////////////////////////////////////////////////////////////////////////
					//add some stuff for the data exchange between the GPUs //////////////////
					if (k == STARTOFFZ)
					{
						parH[level]->sizePlaneSB  += 1;
						if (parH[level]->isSetSendB == false)
						{
							parH[level]->startB = mm[0];
							parH[level]->isSetSendB = true;
						}
					} 
					else if (k == parH[level]->gridNZ + STARTOFFZ - 1)
					{
						parH[level]->sizePlaneST  += 1;
						if (parH[level]->isSetSendT == false)
						{
							parH[level]->startT = mm[0];
							parH[level]->isSetSendT = true;
						}
					}
					else if (k == parH[level]->gridNZ + STARTOFFZ)
					{
						parH[level]->sizePlaneRB  += 1;
						if (parH[level]->isSetRecvB == false)
						{
							parH[level]->endB = mm[0];
							parH[level]->isSetRecvB = true;
						}
					}
					else if (k == STARTOFFZ-1)
					{
						parH[level]->sizePlaneRT  += 1;
						if (parH[level]->isSetRecvT == false)
						{
							parH[level]->endT = mm[0];
							parH[level]->isSetRecvT = true;
						}
					}
					//////////////////////////////////////////////////////////////////////////
					parH[level]->k[mm[0]]    = parH[level]->size_Mat_SP;
					parH[level]->size_Mat_SP = parH[level]->size_Mat_SP + 1;               
					parD[level]->size_Mat_SP = parD[level]->size_Mat_SP + 1;  
				}
				else parH[level]->k[mm[0]] = 0;
			}
		}
	}
	parH[level]->mem_size_real_SP    = sizeof(real     ) * parH[level]->size_Mat_SP;
	parH[level]->mem_size_int_SP        = sizeof(unsigned int) * parH[level]->size_Mat_SP;
	parD[level]->mem_size_real_SP    = sizeof(real     ) * parD[level]->size_Mat_SP;
	parD[level]->mem_size_int_SP        = sizeof(unsigned int) * parD[level]->size_Mat_SP;
}
void Parameter::fillSparse(int level)
{
	unsigned int li = ((parH[level]->gridNX+STARTOFFX-2)-(STARTOFFX+1)-1);
	unsigned int lj = ((parH[level]->gridNY+STARTOFFY-2)-(STARTOFFY+1)-1);
	real globalX, globalY, globalZ;
	//real InitglobalX, InitglobalY, InitglobalZ;
	real PI = 3.141592653589793238462643383279f;

	for (unsigned int k=1; k<parH[level]->gridNZ + 2 * STARTOFFZ - 1; k++)
	{
		for (unsigned int j=1; j<parH[level]->gridNY + 2 * STARTOFFY - 1; j++)
		{
			for (unsigned int i=1; i<parH[level]->gridNX + 2 * STARTOFFX - 1; i++)
			{
				int m = parH[level]->nx*(parH[level]->ny*k + j) + i;
				if ((k < parH[level]->gridNZ + 2 * STARTOFFZ - 2) && (j < parH[level]->gridNY + 2 * STARTOFFY - 2) && (i < parH[level]->gridNX + 2 * STARTOFFX - 2))
				{
					if ((X1PERIODIC == true) && (level==coarse) && (i==parH[level]->gridNX + STARTOFFX - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*k + j) + STARTOFFX;
						parH[level]->neighborX_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborX_SP[parH[level]->k[m]] = parH[level]->k[m+1];
					}
					if ((X2PERIODIC == true) && (level==coarse) && (j==parH[level]->gridNY + STARTOFFY - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*k + STARTOFFY) + i;
						parH[level]->neighborY_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborY_SP[parH[level]->k[m]] = parH[level]->k[m+parH[level]->nx];
					}
					if ((X3PERIODIC == true) && (level==coarse) && (k==parH[level]->gridNZ + STARTOFFZ - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*STARTOFFZ + j) + i;
						parH[level]->neighborZ_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborZ_SP[parH[level]->k[m]] = parH[level]->k[m+(parH[level]->nx*parH[level]->ny)];
					}
				}
				parH[level]->geoSP[parH[level]->k[m]]        = parH[level]->geo[m];
				////////////////////////////////////////////////////////////////////////////
				////Coordinates
				//parH[level]->coordX_SP[parH[level]->k[m]]    = i;
				//parH[level]->coordY_SP[parH[level]->k[m]]    = j;
				//parH[level]->coordZ_SP[parH[level]->k[m]]    = k;
				////////////////////////////////////////////////////////////////////////////
				if (diffOn==true)
				{
					parH[level]->Conc[parH[level]->k[m]]         = parH[level]->Conc_Full[m];
				}
				////////////////////////////////////////////////////////////////////////////
				////set pressure in the middle of the fine grid
				//if (level == getFine())
				//{
				//   if(   i == parH[level]->gridNX/2 + STARTOFFX
				//      && j == parH[level]->gridNY/2 + STARTOFFY 
				//      && k == parH[level]->gridNZ/2 + STARTOFFZ) 
				//      parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.1f;             
				//   else 
				//      parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;
				//} 
				//else
				//{
				//   parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;
				//}
				globalX = TrafoXtoWorld(i,level);
				globalY = TrafoYtoWorld(j,level);
				globalZ = TrafoZtoWorld(k,level);
				//without setting a pressure
				parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;       //parH[level]->Conc_Full[m];//bitte schnell wieder entfernen!!!
				//////////////////////////////////////////////////////////////////////////
				parH[level]->vx_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vx_SP[parH[level]->k[m]]        = u0/3.0;
				parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vy_SP[parH[level]->k[m]]        = u0/3.0;
				parH[level]->vz_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vz_SP[parH[level]->k[m]]        = u0/3.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(u0*2.f)*((-4.f*globalX*globalX + parH[level]->gridNX*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + globalX*(-4.f + 4.f*parH[level]->gridNX + 8.f*STARTOFFX))*(-4.f*globalY*globalY + parH[level]->gridNY*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + globalY*(-4.f + 4.f*parH[level]->gridNY + 8.f*STARTOFFY)))/((2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNY)*(2.f - parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(u0*2.f)*((-4.f*i*i + parH[level]->gridNX*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*parH[level]->gridNX + 8.f*STARTOFFX))*(-4.f*j*j + parH[level]->gridNY*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*parH[level]->gridNY + 8.f*STARTOFFY)))/((2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNY)*(2.f - parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);//(16.f*(u0*2.f)*i*j*(parH[level]->nx-i)*(parH[level]->ny-j))/(parH[level]->nx*parH[level]->nx*parH[level]->ny*parH[level]->ny); //u0;
				//////////////////////////////////////////////////////////////////////////
				////gerade
				//parH[level]->vx_SP[parH[level]->k[m]] = (real)((32. * 32. * 3.) / (1000.*(real)parH[level]->gridNX));//(real)parH[level]->gridNX / (real)1000 * 3.0;
				//parH[level]->vy_SP[parH[level]->k[m]] = (real)((getVelocity() * sin(2.0 * i / parH[level]->gridNX * PI) * cos(2.0 * k / parH[level]->gridNZ * PI)) * (32. / (real)parH[level]->gridNX));
				//parH[level]->vz_SP[parH[level]->k[m]] = (real)0.0f;
				//schräg x
				// 			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + (getVelocity() * cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				// 			parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				// 			parH[level]->vz_SP[parH[level]->k[m]]        = (real)(getVelocity() * cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				//schräg z
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNZ) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));

				  			//Taylor Green Vortex uniform
				  			parH[level]->rho_SP[parH[level]->k[m]]       = (real)((getVelocity()*getVelocity())*3.0/4.0*(cos((i)*4.0*PI/(real)parH[level]->gridNX)+cos((k)*4.0*PI/(real)parH[level]->gridNZ)))*(real)(parH[level]->gridNZ)/(real)(parH[level]->gridNX);
				  			//inkl. überlagerter Geschwindigkeit
				  // 			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000. * 32.) * getVelocity() / 0.001 + getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			//ohne überlagerter Geschwindigkeit
				  //			parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				  			parH[level]->vz_SP[parH[level]->k[m]]        = (real)(-getVelocity()*cos(((i)*2.0*PI/(real)parH[level]->gridNX))*sin((k)*2.0*PI/(real)parH[level]->gridNZ))*(real)(parH[level]->gridNZ)/(real)(parH[level]->gridNX);            

				//Kernel Fix Test
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				////parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				////parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				////parH[level]->vz_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNZ) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				//////////////////////////////////////////////////////////////////////////
				//Taylor Green Vortex
				//InitglobalX = TrafoXtoMGsWorld(i,level);
				//InitglobalY = TrafoYtoMGsWorld(j,level);
				//InitglobalZ = TrafoZtoMGsWorld(k,level);
				//parH[level]->rho_SP[parH[level]->k[m]]       = (real)((u0*u0)*3.f/4.f*(cos((InitglobalX)*4.f*PI/parH[level]->gridNX)+cos((InitglobalY)*4.f*PI/parH[level]->gridNY)));
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)( u0*sin(((InitglobalX)*2.f*PI/parH[level]->gridNX))*cos((InitglobalY)*2.f*PI/parH[level]->gridNY));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)(-u0*cos(((InitglobalX)*2.f*PI/parH[level]->gridNX))*sin((InitglobalY)*2.f*PI/parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)0.0f;            
				//////////////////////////////////////////////////////////////////////////
			}
		}
	}
	parH[level]->neighborX_SP[parH[level]->k[0]] = 0;
	parH[level]->neighborY_SP[parH[level]->k[0]] = 0;
	parH[level]->neighborZ_SP[parH[level]->k[0]] = 0;
	parH[level]->geoSP[       parH[level]->k[0]] = GEO_VOID;
	parH[level]->rho_SP[      parH[level]->k[0]] = (real)0.f;
	parH[level]->vx_SP[       parH[level]->k[0]] = (real)0.f;
	parH[level]->vy_SP[       parH[level]->k[0]] = (real)0.f;
	parH[level]->vz_SP[       parH[level]->k[0]] = (real)0.f;
	////////////////////////////////////////////////////////////////////////////
	////Coordinates
	//parH[level]->coordX_SP[parH[level]->k[0]]    = 0;
	//parH[level]->coordY_SP[parH[level]->k[0]]    = 0;
	//parH[level]->coordZ_SP[parH[level]->k[0]]    = 0;
	////////////////////////////////////////////////////////////////////////////
}
void Parameter::copyMeasurePointsArrayToVector(int lev)
{
	int valuesPerClockCycle = (int)(getclockCycleForMP()/getTimestepForMP());
	for(int i = 0; i < (int)parH[lev]->MP.size(); i++)
	{
		for(int j = 0; j < valuesPerClockCycle; j++)
		{
			int index = i*valuesPerClockCycle+j;
			parH[lev]->MP[i].Vx.push_back(parH[lev]->VxMP[index]);
			parH[lev]->MP[i].Vy.push_back(parH[lev]->VyMP[index]);
			parH[lev]->MP[i].Vz.push_back(parH[lev]->VzMP[index]);
			parH[lev]->MP[i].Rho.push_back(parH[lev]->RhoMP[index]);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cuda-alloc-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//full
void Parameter::cudaAllocFull(int lev)
{
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->geo      ), parH[lev]->mem_size_int  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->k        ), parH[lev]->mem_size_int  ));
}
void Parameter::cudaFreeFull(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->geo   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->k     ));
}
//coord
void Parameter::cudaAllocCoord(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->coordX_SP      ), parH[lev]->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->coordY_SP      ), parH[lev]->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->coordZ_SP      ), parH[lev]->mem_size_real_SP  ));
	//Device (spinning ship)
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->coordX_SP      ), parH[lev]->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->coordY_SP      ), parH[lev]->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->coordZ_SP      ), parH[lev]->mem_size_real_SP  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
	//printf("Coord = %f MB",tmp/1000000.);  
}
void Parameter::cudaCopyCoord(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parD[lev]->coordX_SP,  parH[lev]->coordX_SP,  parH[lev]->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->coordY_SP,  parH[lev]->coordY_SP,  parH[lev]->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->coordZ_SP,  parH[lev]->coordZ_SP,  parH[lev]->mem_size_real_SP     , cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeCoord(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->coordX_SP   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->coordY_SP   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->coordZ_SP   ));
}
//print
void Parameter::cudaCopyPrint(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->vx_SP   , parD[lev]->vx_SP   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->vy_SP   , parD[lev]->vy_SP   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->vz_SP   , parD[lev]->vz_SP   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->rho_SP  , parD[lev]->rho_SP  , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->press_SP, parD[lev]->press_SP, parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void Parameter::cudaCopyMedianPrint(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->vx_SP_Med   , parD[lev]->vx_SP_Med   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->vy_SP_Med   , parD[lev]->vy_SP_Med   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->vz_SP_Med   , parD[lev]->vz_SP_Med   , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->rho_SP_Med  , parD[lev]->rho_SP_Med  , parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->press_SP_Med, parD[lev]->press_SP_Med, parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
//sparse
void Parameter::cudaAllocSP(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->geoSP           ), parH[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->neighborX_SP    ), parH[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->neighborY_SP    ), parH[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->neighborZ_SP    ), parH[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->rho_SP          ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vx_SP           ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vy_SP           ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vz_SP           ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->press_SP        ), parH[lev]->mem_size_real_SP));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->geoSP               ), parD[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->neighborX_SP        ), parD[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->neighborY_SP        ), parD[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->neighborZ_SP        ), parD[lev]->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->rho_SP              ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vx_SP               ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vy_SP               ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vz_SP               ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->press_SP            ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->d0SP.f[0]           ), (unsigned long long)getD3Qxx()*(unsigned long long)parD[lev]->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	//double tmp = 4. * (double)parH[lev]->mem_size_int_SP + 10. * (double)parH[lev]->mem_size_real_SP + (double)getD3Qxx() * (double)parH[lev]->mem_size_real_SP;
	double tmp = 4. * (double)parH[lev]->mem_size_int_SP + 5. * (double)parH[lev]->mem_size_real_SP + (double)getD3Qxx() * (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);

	//int test = sizeof(int);
	//printf("\n sizeof int = %d \n",test); 

	//printf("AlocSP = %f MB \n",tmp/1000000.);  
	//int test = sizeof(float*);
	//printf("float* = %d \n",test); 
	//unsigned long long test2 = (unsigned long long)getD3Qxx()*(unsigned long long)parD[lev]->mem_size_real_SP;
	//test2 = test2 / 1000000000; 
	//printf("test2 = %d \n",test2); 
}
void Parameter::cudaCopySP(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parD[lev]->geoSP       ,  parH[lev]->geoSP       ,  parH[lev]->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->neighborX_SP,  parH[lev]->neighborX_SP,  parH[lev]->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->neighborY_SP,  parH[lev]->neighborY_SP,  parH[lev]->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->neighborZ_SP,  parH[lev]->neighborZ_SP,  parH[lev]->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->rho_SP      ,  parH[lev]->rho_SP      ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vx_SP       ,  parH[lev]->vx_SP       ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vy_SP       ,  parH[lev]->vy_SP       ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vz_SP       ,  parH[lev]->vz_SP       ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->press_SP    ,  parH[lev]->press_SP    ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeSP(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->geoSP       ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vx_SP       ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vy_SP       ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vz_SP       ));
	checkCudaErrors( cudaFreeHost(parH[lev]->rho_SP      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->press_SP    ));
	checkCudaErrors( cudaFreeHost(parH[lev]->neighborX_SP));
	checkCudaErrors( cudaFreeHost(parH[lev]->neighborY_SP));
	checkCudaErrors( cudaFreeHost(parH[lev]->neighborZ_SP));
}
//F3
void Parameter::cudaAllocF3SP(int lev)
{
	//Device						 
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->g6.g[0]), (unsigned long long)6*(unsigned long long)parD[lev]->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)6 * (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
//negative neighbor (WSB)
void Parameter::cudaAllocNeighborWSB(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->neighborWSB_SP    ), parH[lev]->mem_size_int_SP    ));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->neighborWSB_SP        ), parD[lev]->mem_size_int_SP    ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->mem_size_int_SP;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyNeighborWSB(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parD[lev]->neighborWSB_SP,  parH[lev]->neighborWSB_SP,  parH[lev]->mem_size_int_SP     , cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeNeighborWSB(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->neighborWSB_SP));
}
//turbulent viscosity
void Parameter::cudaAllocTurbulentViscosity(int lev)
{
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->turbViscosity), parH[lev]->mem_size_real_SP));
	//Debug
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gSij ), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gSDij), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDxvx), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDyvx), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDzvx), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDxvy), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDyvy), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDzvy), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDxvz), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDyvz), parH[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parH[lev]->gDzvz), parH[lev]->mem_size_real_SP));

	//Device						 
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->turbViscosity), parD[lev]->mem_size_real_SP));
	//Debug
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gSij ), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gSDij), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDxvx), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDyvx), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDzvx), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDxvy), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDyvy), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDzvy), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDxvz), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDyvz), parD[lev]->mem_size_real_SP));
	checkCudaErrors(cudaMalloc((void**) &(parD[lev]->gDzvz), parD[lev]->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyTurbulentViscosityHD(int lev)
{
	//copy host to device
	checkCudaErrors(cudaMemcpy(parD[lev]->turbViscosity, parH[lev]->turbViscosity, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	//Debug
	checkCudaErrors(cudaMemcpy(parD[lev]->gSij , parH[lev]->gSij , parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gSDij, parH[lev]->gSDij, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDxvx, parH[lev]->gDxvx, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDyvx, parH[lev]->gDyvx, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDzvx, parH[lev]->gDzvx, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDxvy, parH[lev]->gDxvy, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDyvy, parH[lev]->gDyvy, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDzvy, parH[lev]->gDzvy, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDxvz, parH[lev]->gDxvz, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDyvz, parH[lev]->gDyvz, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(parD[lev]->gDzvz, parH[lev]->gDzvz, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyTurbulentViscosityDH(int lev)
{
	//copy device to host
	checkCudaErrors(cudaMemcpy(parH[lev]->turbViscosity, parD[lev]->turbViscosity, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	//Debug
	checkCudaErrors(cudaMemcpy(parH[lev]->gSij , parD[lev]->gSij , parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gSDij, parD[lev]->gSDij, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDxvx, parD[lev]->gDxvx, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDyvx, parD[lev]->gDyvx, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDzvx, parD[lev]->gDzvx, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDxvy, parD[lev]->gDxvy, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDyvy, parD[lev]->gDyvy, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDzvy, parD[lev]->gDzvy, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDxvz, parD[lev]->gDxvz, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDyvz, parD[lev]->gDyvz, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(parH[lev]->gDzvz, parD[lev]->gDzvz, parH[lev]->mem_size_real_SP, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeTurbulentViscosity(int lev)
{
	checkCudaErrors(cudaFreeHost(parH[lev]->turbViscosity));
	//Debug
	checkCudaErrors(cudaFreeHost(parH[lev]->gSij ));
	checkCudaErrors(cudaFreeHost(parH[lev]->gSDij));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDxvx));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDyvx));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDzvx));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDxvy));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDyvy));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDzvy));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDxvz));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDyvz));
	checkCudaErrors(cudaFreeHost(parH[lev]->gDzvz));
}
//median
void Parameter::cudaAllocMedianSP(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->rho_SP_Med      ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vx_SP_Med       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vy_SP_Med       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vz_SP_Med       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->press_SP_Med    ), parH[lev]->mem_size_real_SP));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->rho_SP_Med          ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vx_SP_Med           ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vy_SP_Med           ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->vz_SP_Med           ), parD[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->press_SP_Med        ), parD[lev]->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 5. * (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyMedianSP(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parD[lev]->rho_SP_Med  ,  parH[lev]->rho_SP_Med  ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vx_SP_Med   ,  parH[lev]->vx_SP_Med   ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vy_SP_Med   ,  parH[lev]->vy_SP_Med   ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->vz_SP_Med   ,  parH[lev]->vz_SP_Med   ,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->press_SP_Med,  parH[lev]->press_SP_Med,  parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeMedianSP(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->vx_SP_Med   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vy_SP_Med   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vz_SP_Med   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->rho_SP_Med  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->press_SP_Med));
}
void Parameter::cudaAllocMedianOut(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->rho_SP_Med_Out      ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vx_SP_Med_Out       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vy_SP_Med_Out       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->vz_SP_Med_Out       ), parH[lev]->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->press_SP_Med_Out    ), parH[lev]->mem_size_real_SP));
}
void Parameter::cudaFreeMedianOut(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->vx_SP_Med_Out   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vy_SP_Med_Out   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->vz_SP_Med_Out   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->rho_SP_Med_Out  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->press_SP_Med_Out));
}
//Interface CF
void Parameter::cudaAllocInterfaceCF(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->intCF.ICellCFC), parH[lev]->mem_size_kCF  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->intCF.ICellCFF), parH[lev]->mem_size_kCF  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->intCF.ICellCFC), parD[lev]->mem_size_kCF  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->intCF.ICellCFF), parD[lev]->mem_size_kCF  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)parH[lev]->mem_size_kCF;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInterfaceCF(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->intCF.ICellCFC, parH[lev]->intCF.ICellCFC, parH[lev]->mem_size_kCF, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->intCF.ICellCFF, parH[lev]->intCF.ICellCFF, parH[lev]->mem_size_kCF, cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeInterfaceCF(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->intCF.ICellCFC));
	checkCudaErrors( cudaFreeHost(parH[lev]->intCF.ICellCFF));
}
//Interface FC
void Parameter::cudaAllocInterfaceFC(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->intFC.ICellFCF), parH[lev]->mem_size_kFC  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->intFC.ICellFCC), parH[lev]->mem_size_kFC  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->intFC.ICellFCF), parD[lev]->mem_size_kFC  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->intFC.ICellFCC), parD[lev]->mem_size_kFC  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)parH[lev]->mem_size_kFC;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInterfaceFC(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->intFC.ICellFCF, parH[lev]->intFC.ICellFCF, parH[lev]->mem_size_kFC, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->intFC.ICellFCC, parH[lev]->intFC.ICellFCC, parH[lev]->mem_size_kFC, cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeInterfaceFC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->intFC.ICellFCF));
	checkCudaErrors( cudaFreeHost(parH[lev]->intFC.ICellFCC));
}
//Interface Offset CF
void Parameter::cudaAllocInterfaceOffCF(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offCF.xOffCF),   parH[lev]->mem_size_kCF_off  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offCF.yOffCF),   parH[lev]->mem_size_kCF_off  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offCF.zOffCF),   parH[lev]->mem_size_kCF_off  ));
	getLastCudaError("Allocate host memory");
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offCF.xOffCF),   parD[lev]->mem_size_kCF_off  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offCF.yOffCF),   parD[lev]->mem_size_kCF_off  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offCF.zOffCF),   parD[lev]->mem_size_kCF_off  ));
	getLastCudaError("Allocate device memory");
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parH[lev]->mem_size_kCF_off;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInterfaceOffCF(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->offCF.xOffCF,   parH[lev]->offCF.xOffCF,   parH[lev]->mem_size_kCF_off, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->offCF.yOffCF,   parH[lev]->offCF.yOffCF,   parH[lev]->mem_size_kCF_off, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->offCF.zOffCF,   parH[lev]->offCF.zOffCF,   parH[lev]->mem_size_kCF_off, cudaMemcpyHostToDevice));
	getLastCudaError("Copy host memory to device");
}
void Parameter::cudaFreeInterfaceOffCF(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->offCF.xOffCF));
	checkCudaErrors( cudaFreeHost(parH[lev]->offCF.yOffCF));
	checkCudaErrors( cudaFreeHost(parH[lev]->offCF.zOffCF));
}
//Interface Offset FC
void Parameter::cudaAllocInterfaceOffFC(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offFC.xOffFC),   parH[lev]->mem_size_kFC_off  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offFC.yOffFC),   parH[lev]->mem_size_kFC_off  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->offFC.zOffFC),   parH[lev]->mem_size_kFC_off  ));
	getLastCudaError("Allocate host memory");
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offFC.xOffFC),   parD[lev]->mem_size_kFC_off  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offFC.yOffFC),   parD[lev]->mem_size_kFC_off  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->offFC.zOffFC),   parD[lev]->mem_size_kFC_off  ));
	getLastCudaError("Allocate device memory");
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parH[lev]->mem_size_kFC_off;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInterfaceOffFC(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->offFC.xOffFC,   parH[lev]->offFC.xOffFC,   parH[lev]->mem_size_kFC_off, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->offFC.yOffFC,   parH[lev]->offFC.yOffFC,   parH[lev]->mem_size_kFC_off, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->offFC.zOffFC,   parH[lev]->offFC.zOffFC,   parH[lev]->mem_size_kFC_off, cudaMemcpyHostToDevice));
	getLastCudaError("Copy host memory to device");
}
void Parameter::cudaFreeInterfaceOffFC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->offFC.xOffFC));
	checkCudaErrors( cudaFreeHost(parH[lev]->offFC.yOffFC));
	checkCudaErrors( cudaFreeHost(parH[lev]->offFC.zOffFC));
}

//Velo
void Parameter::cudaAllocVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parH[lev]->Qinflow.kArray;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parH[lev]->Qinflow.kArray;
	//unsigned int mem_size_inflow_Q_k = sizeof(int)*parH[lev]->Qinflow.kQ;
	//unsigned int mem_size_inflow_Q_q = sizeof(real)*parH[lev]->Qinflow.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.q27[0]),  getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.k),                  mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.Vx),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.Vy),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.Vz),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.deltaVz),            mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qinflow.RhoBC),              mem_size_inflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.q27[0]),      getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.k),                      mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.Vx),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.Vy),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.Vz),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qinflow.deltaVz),                mem_size_inflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_inflow_Q_k + 4. * (double)mem_size_inflow_Q_q + (double)getD3Qxx() * (double)mem_size_inflow_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parH[lev]->Qinflow.kArray;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parH[lev]->Qinflow.kArray;
	//unsigned int mem_size_inflow_Q_k = sizeof(int)*parH[lev]->Qinflow.kQ;
	//unsigned int mem_size_inflow_Q_q = sizeof(real)*parH[lev]->Qinflow.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.q27[0],  parH[lev]->Qinflow.q27[0], getD3Qxx()* mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.k,       parH[lev]->Qinflow.k,                  mem_size_inflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.Vx,      parH[lev]->Qinflow.Vx,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.Vy,      parH[lev]->Qinflow.Vy,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.Vz,      parH[lev]->Qinflow.Vz,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qinflow.deltaVz, parH[lev]->Qinflow.deltaVz,            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));

}
void Parameter::cudaFreeVeloBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.Vx     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.Vy     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.Vz     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qinflow.deltaVz));
}
//Press
void Parameter::cudaAllocOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parH[lev]->Qoutflow.kQ;
	unsigned int mem_size_outflow_Q_q = sizeof(real)*parH[lev]->Qoutflow.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qoutflow.q27[0]), getD3Qxx()*mem_size_outflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qoutflow.k),                 mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qoutflow.kN),                mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Qoutflow.RhoBC),             mem_size_outflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qoutflow.q27[0]),     getD3Qxx()* mem_size_outflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qoutflow.k),                      mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qoutflow.kN),                     mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Qoutflow.RhoBC),                  mem_size_outflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_outflow_Q_k + 2. * (double)mem_size_outflow_Q_q + (double)getD3Qxx()*(double)mem_size_outflow_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parH[lev]->Qoutflow.kQ;
	unsigned int mem_size_outflow_Q_q = sizeof(real)*parH[lev]->Qoutflow.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->Qoutflow.q27[0],  parH[lev]->Qoutflow.q27[0], getD3Qxx()* mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qoutflow.k,       parH[lev]->Qoutflow.k,                  mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qoutflow.kN,      parH[lev]->Qoutflow.kN,                 mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Qoutflow.RhoBC,   parH[lev]->Qoutflow.RhoBC,              mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeOutflowBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->Qoutflow.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qoutflow.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qoutflow.kN     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->Qoutflow.RhoBC  ));
}
//Inlet
void Parameter::cudaAllocInlet(int lev)
{
	unsigned int mem_size_inlet_Q_k = sizeof(int)*parH[lev]->QInlet.kQ;
	unsigned int mem_size_inlet_Q_q = sizeof(real)*parH[lev]->QInlet.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInlet.q27[0]), getD3Qxx()*mem_size_inlet_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInlet.k),                 mem_size_inlet_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInlet.kN),                mem_size_inlet_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInlet.RhoBC),             mem_size_inlet_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInlet.q27[0]),     getD3Qxx()* mem_size_inlet_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInlet.k),                      mem_size_inlet_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInlet.kN),                     mem_size_inlet_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInlet.RhoBC),                  mem_size_inlet_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_inlet_Q_k + (double)mem_size_inlet_Q_q + (double)getD3Qxx()*(double)mem_size_inlet_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInlet(int lev)
{
	unsigned int mem_size_inlet_Q_k = sizeof(int)*parH[lev]->QInlet.kQ;
	unsigned int mem_size_inlet_Q_q = sizeof(real)*parH[lev]->QInlet.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QInlet.q27[0],  parH[lev]->QInlet.q27[0], getD3Qxx()* mem_size_inlet_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInlet.k,       parH[lev]->QInlet.k,                  mem_size_inlet_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInlet.kN,      parH[lev]->QInlet.kN,                 mem_size_inlet_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInlet.RhoBC,   parH[lev]->QInlet.RhoBC,              mem_size_inlet_Q_q,  cudaMemcpyHostToDevice));
}																  
void Parameter::cudaFreeInlet(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QInlet.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInlet.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInlet.kN     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInlet.RhoBC  ));
}
//Outlet
void Parameter::cudaAllocOutlet(int lev)
{
	unsigned int mem_size_outlet_Q_k = sizeof(int)*parH[lev]->QOutlet.kQ;
	unsigned int mem_size_outlet_Q_q = sizeof(real)*parH[lev]->QOutlet.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutlet.q27[0]), getD3Qxx()*mem_size_outlet_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutlet.k),                 mem_size_outlet_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutlet.kN),                mem_size_outlet_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutlet.RhoBC),             mem_size_outlet_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutlet.q27[0]),     getD3Qxx()* mem_size_outlet_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutlet.k),                      mem_size_outlet_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutlet.kN),                     mem_size_outlet_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutlet.RhoBC),                  mem_size_outlet_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_outlet_Q_k + (double)mem_size_outlet_Q_q + (double)getD3Qxx()*(double)mem_size_outlet_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyOutlet(int lev)
{
	unsigned int mem_size_outlet_Q_k = sizeof(int)*parH[lev]->QOutlet.kQ;
	unsigned int mem_size_outlet_Q_q = sizeof(real)*parH[lev]->QOutlet.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QOutlet.q27[0],  parH[lev]->QOutlet.q27[0], getD3Qxx()* mem_size_outlet_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutlet.k,       parH[lev]->QOutlet.k,                  mem_size_outlet_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutlet.kN,      parH[lev]->QOutlet.kN,                 mem_size_outlet_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutlet.RhoBC,   parH[lev]->QOutlet.RhoBC,              mem_size_outlet_Q_q,  cudaMemcpyHostToDevice));
}																  
void Parameter::cudaFreeOutlet(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutlet.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutlet.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutlet.kN     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutlet.RhoBC  ));
}
//Wall
void Parameter::cudaAllocWallBC(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QWall.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QWall.kQ;
	unsigned int mem_size_Q_value  = sizeof(long long)*parH[lev]->QWall.kQ; //Geller
	unsigned int mem_size_Q_q_read = sizeof(real)*parH[lev]->kQread;     //Geller

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QWall.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QWall.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QWall.qread),             mem_size_Q_q_read ));//Geller
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QWall.valueQ),            mem_size_Q_value  ));//Geller

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QWall.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QWall.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyWallBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QWall.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QWall.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QWall.q27[0], parH[lev]->QWall.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QWall.k,      parH[lev]->QWall.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeWallBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QWall.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QWall.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QWall.valueQ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QWall.qread));
}
//Geometrie
void Parameter::cudaAllocGeomBC(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QGeom.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QGeom.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeom.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeom.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeom.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeom.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyGeomBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QGeom.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QGeom.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QGeom.q27[0], parH[lev]->QGeom.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeom.k,      parH[lev]->QGeom.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeGeomBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeom.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeom.k));
}
//Geometrie inkl. Values
void Parameter::cudaAllocGeomValuesBC(int lev)
{
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QGeom.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeom.Vx),  mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeom.Vy),  mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeom.Vz),  mem_size_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeom.Vx),      mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeom.Vy),      mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeom.Vz),      mem_size_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyGeomValuesBC(int lev)
{
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QGeom.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QGeom.Vx, parH[lev]->QGeom.Vx,  mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeom.Vy, parH[lev]->QGeom.Vy,  mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeom.Vz, parH[lev]->QGeom.Vz,  mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeGeomValuesBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeom.Vx));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeom.Vy));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeom.Vz));
}
//Geometrie inkl. Normale für Slip
void Parameter::cudaAllocGeomNormals(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QGeomNormalX.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QGeomNormalX.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalX.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalX.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalY.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalY.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalZ.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QGeomNormalZ.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalX.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalX.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalY.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalY.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalZ.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QGeomNormalZ.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyGeomNormals(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QGeomNormalX.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QGeomNormalX.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalX.q27[0], parH[lev]->QGeomNormalX.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalX.k,      parH[lev]->QGeomNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalY.q27[0], parH[lev]->QGeomNormalY.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalY.k,      parH[lev]->QGeomNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalZ.q27[0], parH[lev]->QGeomNormalZ.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QGeomNormalZ.k,      parH[lev]->QGeomNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeGeomNormals(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalX.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalX.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalY.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalY.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalZ.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QGeomNormalZ.k));
}
//Geometrie inkl. Normale für Inflow
void Parameter::cudaAllocInflowNormals(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QInflowNormalX.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QInflowNormalX.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalX.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalX.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalY.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalY.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalZ.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QInflowNormalZ.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalX.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalX.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalY.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalY.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalZ.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QInflowNormalZ.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyInflowNormals(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QInflowNormalX.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QInflowNormalX.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalX.q27[0], parH[lev]->QInflowNormalX.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalX.k,      parH[lev]->QInflowNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalY.q27[0], parH[lev]->QInflowNormalY.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalY.k,      parH[lev]->QInflowNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalZ.q27[0], parH[lev]->QInflowNormalZ.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QInflowNormalZ.k,      parH[lev]->QInflowNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeInflowNormals(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalX.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalX.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalY.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalY.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalZ.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QInflowNormalZ.k));
}
//Geometrie inkl. Normale für Outflow
void Parameter::cudaAllocOutflowNormals(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QOutflowNormalX.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QOutflowNormalX.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalX.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalX.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalY.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalY.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalZ.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QOutflowNormalZ.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalX.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalX.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalY.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalY.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalZ.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QOutflowNormalZ.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyOutflowNormals(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QOutflowNormalX.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QOutflowNormalX.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalX.q27[0], parH[lev]->QOutflowNormalX.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalX.k,      parH[lev]->QOutflowNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalY.q27[0], parH[lev]->QOutflowNormalY.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalY.k,      parH[lev]->QOutflowNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalZ.q27[0], parH[lev]->QOutflowNormalZ.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QOutflowNormalZ.k,      parH[lev]->QOutflowNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeOutflowNormals(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalX.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalX.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalY.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalY.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalZ.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QOutflowNormalZ.k));
}
//Slip
void Parameter::cudaAllocSlipBC(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QSlip.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QSlip.kQ;
	//unsigned int mem_size_Q_value  = sizeof(long long)*parH[lev]->QSlip.kQ; //Geller
	//unsigned int mem_size_Q_q_read = sizeof(real)*parH[lev]->kSlipQread;     //Geller

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QSlip.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QSlip.k),                 mem_size_Q_k      ));
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QSlip.qread),             mem_size_Q_q_read ));//Geller
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QSlip.valueQ),            mem_size_Q_value  ));//Geller

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QSlip.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QSlip.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopySlipBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QSlip.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QSlip.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QSlip.q27[0], parH[lev]->QSlip.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QSlip.k,      parH[lev]->QSlip.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeSlipBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QSlip.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QSlip.k));
	//checkCudaErrors( cudaFreeHost(parH[lev]->QSlip.valueQ));
	//checkCudaErrors( cudaFreeHost(parH[lev]->QSlip.qread));
}
//Press (no Geller)
void Parameter::cudaAllocPress(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parH[lev]->QPress.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parH[lev]->QPress.kQ;
	//unsigned int mem_size_Q_value  = sizeof(long long)*parH[lev]->QPress.kQ; //Geller
	//unsigned int mem_size_Q_q_read = sizeof(real)*parH[lev]->kPressQread;     //Geller

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.q27[0]), getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.kN),                mem_size_Q_k      )); 
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.RhoBC),             mem_size_Q_q      ));
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.qread),             mem_size_Q_q_read ));//Geller
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPress.valueQ),            mem_size_Q_value  ));//Geller

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPress.q27[0]),     getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPress.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPress.kN),                     mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPress.RhoBC),                  mem_size_Q_q     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_Q_k + (double)mem_size_Q_q + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPress(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QPress.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QPress.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QPress.q27[0], parH[lev]->QPress.q27[0], getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPress.k,      parH[lev]->QPress.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPress.kN,     parH[lev]->QPress.kN,                 mem_size_Q_k,       cudaMemcpyHostToDevice)); 
	checkCudaErrors( cudaMemcpy(parD[lev]->QPress.RhoBC,  parH[lev]->QPress.RhoBC,              mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void Parameter::cudaFreePress(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QPress.q27[0]));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPress.k));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPress.kN));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPress.RhoBC));
	//checkCudaErrors( cudaFreeHost(parH[lev]->QPress.valueQ));//Geller
	//checkCudaErrors( cudaFreeHost(parH[lev]->QPress.qread));//Geller
}
//Test roundoff error
void Parameter::cudaAllocTestRE(int lev, unsigned int size)
{
	unsigned int mem_size = sizeof(real)*size;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kDistTestRE.f[0]), (1+getD3Qxx())*mem_size));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kDistTestRE.f[0]), (1+getD3Qxx())*mem_size));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyTestREtoDevice(int lev, unsigned int size)
{
	unsigned int mem_size = sizeof(real)*size;
	checkCudaErrors( cudaMemcpy(parD[lev]->kDistTestRE.f[0], parH[lev]->kDistTestRE.f[0], (1+getD3Qxx())*mem_size, cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyTestREtoHost(int lev, unsigned int size)
{
	unsigned int mem_size = sizeof(real)*size;
	checkCudaErrors( cudaMemcpy(parH[lev]->kDistTestRE.f[0], parD[lev]->kDistTestRE.f[0], (1+getD3Qxx())*mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeTestRE(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->kDistTestRE.f[0]));
}
//PressX0 = X-inflow
void Parameter::cudaAllocPressX0(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QpressX0.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QpressX0.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.q27[0]), getD3Qxx()*mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.k),                 mem_size_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.kN),                mem_size_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.Vx),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.Vy),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.Vz),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.RhoBC),             mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX0.deltaVz),           mem_size_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.q27[0]),     getD3Qxx()* mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.k),                      mem_size_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.kN),                     mem_size_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.Vx),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.Vy),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.Vz),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.RhoBC),                  mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX0.deltaVz),                mem_size_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_Q_k + 5. * (double)mem_size_Q_q + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPressX0(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QpressX0.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QpressX0.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.q27[0],  parH[lev]->QpressX0.q27[0], getD3Qxx()* mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.k,       parH[lev]->QpressX0.k,                  mem_size_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.kN,      parH[lev]->QpressX0.kN,                 mem_size_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.Vx,      parH[lev]->QpressX0.Vx,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.Vy,      parH[lev]->QpressX0.Vy,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.Vz,      parH[lev]->QpressX0.Vz,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.RhoBC,   parH[lev]->QpressX0.RhoBC,              mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX0.deltaVz, parH[lev]->QpressX0.deltaVz,            mem_size_Q_q,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreePressX0(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.kN     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.Vx     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.Vy     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.Vz     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.RhoBC  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX0.deltaVz));
}
//PressX1 = X-outflow
void Parameter::cudaAllocPressX1(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QpressX1.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QpressX1.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.q27[0]), getD3Qxx()*mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.k),                 mem_size_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.kN),                mem_size_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.Vx),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.Vy),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.Vz),                mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.RhoBC),             mem_size_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QpressX1.deltaVz),           mem_size_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.q27[0]),     getD3Qxx()* mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.k),                      mem_size_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.kN),                     mem_size_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.Vx),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.Vy),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.Vz),                     mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.RhoBC),                  mem_size_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QpressX1.deltaVz),                mem_size_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_Q_k + 5. * (double)mem_size_Q_q + (double)getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPressX1(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parH[lev]->QpressX1.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parH[lev]->QpressX1.kQ;

	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.q27[0],  parH[lev]->QpressX1.q27[0], getD3Qxx()* mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.k,       parH[lev]->QpressX1.k,                  mem_size_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.kN,      parH[lev]->QpressX1.kN,                 mem_size_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.Vx,      parH[lev]->QpressX1.Vx,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.Vy,      parH[lev]->QpressX1.Vy,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.Vz,      parH[lev]->QpressX1.Vz,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.RhoBC,   parH[lev]->QpressX1.RhoBC,              mem_size_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QpressX1.deltaVz, parH[lev]->QpressX1.deltaVz,            mem_size_Q_q,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreePressX1(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.kN     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.Vx     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.Vy     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.Vz     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.RhoBC  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QpressX1.deltaVz));
}
//Propeller Velocity
void Parameter::cudaAllocVeloPropeller(int lev)
{
	unsigned int mem_size_Propeller_k = sizeof(int)*parH[lev]->QPropeller.kQ;
	unsigned int mem_size_Propeller_q = sizeof(real)*parH[lev]->QPropeller.kQ;

	//Host
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.q27[0]),  getD3Qxx()*mem_size_Propeller_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.k),                  mem_size_Propeller_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.Vx),                 mem_size_Propeller_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.Vy),                 mem_size_Propeller_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.Vz),                 mem_size_Propeller_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->QPropeller.RhoBC),              mem_size_Propeller_q ));

	//Device
	//checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.q27[0]),      getD3Qxx()*mem_size_Propeller_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.k),                      mem_size_Propeller_k ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.Vx),                     mem_size_Propeller_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.Vy),                     mem_size_Propeller_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.Vz),                     mem_size_Propeller_q ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->QPropeller.RhoBC),                  mem_size_Propeller_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Propeller_k + 4. * (double)mem_size_Propeller_q;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyVeloPropeller(int lev)
{
	unsigned int mem_size_Propeller_k = sizeof(int)*parH[lev]->QPropeller.kQ;
	unsigned int mem_size_Propeller_q = sizeof(real)*parH[lev]->QPropeller.kQ;

	//checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.q27[0],  parH[lev]->QPropeller.q27[0], getD3Qxx()* mem_size_Propeller_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.k,       parH[lev]->QPropeller.k,                  mem_size_Propeller_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.Vx,      parH[lev]->QPropeller.Vx,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.Vy,      parH[lev]->QPropeller.Vy,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.Vz,      parH[lev]->QPropeller.Vz,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->QPropeller.RhoBC,   parH[lev]->QPropeller.RhoBC,              mem_size_Propeller_q,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeVeloPropeller(int lev)
{
	//checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.q27[0] ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.k      ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.Vx     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.Vy     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.Vz     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->QPropeller.RhoBC  ));
}
//Measure Points
//void Parameter::cudaAllocMeasurePoints(int lev, int i)
//{
//	//Host
//	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->MP[i].Vx),                 parH[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->MP[i].Vy),                 parH[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->MP[i].Vz),                 parH[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->MP[i].Rho),                parH[lev]->memSizerealMP ));
//
//	//Device
//	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->MP[i].Vx),                     parD[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->MP[i].Vy),                     parD[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->MP[i].Vz),                     parD[lev]->memSizerealMP ));
//	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->MP[i].Rho),                    parD[lev]->memSizerealMP ));
//}
//void Parameter::cudaCopyMeasurePoints(int lev, int i)
//{
//	checkCudaErrors( cudaMemcpy(parD[lev]->MP[i].Vx,      parH[lev]->MP[i].Vx,           parH[lev]->memSizerealMP,  cudaMemcpyHostToDevice));
//	checkCudaErrors( cudaMemcpy(parD[lev]->MP[i].Vy,      parH[lev]->MP[i].Vy,           parH[lev]->memSizerealMP,  cudaMemcpyHostToDevice));
//	checkCudaErrors( cudaMemcpy(parD[lev]->MP[i].Vz,      parH[lev]->MP[i].Vz,           parH[lev]->memSizerealMP,  cudaMemcpyHostToDevice));
//	checkCudaErrors( cudaMemcpy(parD[lev]->MP[i].Rho,     parH[lev]->MP[i].Rho,          parH[lev]->memSizerealMP,  cudaMemcpyHostToDevice));
//}
//void Parameter::cudaFreeMeasurePoints(int lev, int i)
//{
//	checkCudaErrors( cudaFreeHost(parH[lev]->MP[i].Vx     ));
//	checkCudaErrors( cudaFreeHost(parH[lev]->MP[i].Vy     ));
//	checkCudaErrors( cudaFreeHost(parH[lev]->MP[i].Vz     ));
//	checkCudaErrors( cudaFreeHost(parH[lev]->MP[i].Rho    ));
//}
void Parameter::cudaAllocMeasurePointsIndex(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kMP),						parH[lev]->memSizeIntkMP     ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VxMP),					parH[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VyMP),					parH[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VzMP),					parH[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->RhoMP),					parH[lev]->memSizerealkMP ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kMP),							parD[lev]->memSizeIntkMP     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VxMP),						parD[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VyMP),						parD[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VzMP),						parD[lev]->memSizerealkMP ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->RhoMP),						parD[lev]->memSizerealkMP ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->memSizeIntkMP + 4. * (double)parH[lev]->memSizerealkMP;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyMeasurePointsIndex(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->kMP,           parH[lev]->kMP,                parH[lev]->memSizeIntkMP,      cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->VxMP,          parH[lev]->VxMP,               parH[lev]->memSizerealkMP,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->VyMP,          parH[lev]->VyMP,               parH[lev]->memSizerealkMP,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->VzMP,          parH[lev]->VzMP,               parH[lev]->memSizerealkMP,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->RhoMP,         parH[lev]->RhoMP,              parH[lev]->memSizerealkMP,  cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyMeasurePointsToHost(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->kMP,           parD[lev]->kMP,                parH[lev]->memSizeIntkMP,      cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->VxMP,          parD[lev]->VxMP,               parH[lev]->memSizerealkMP,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->VyMP,          parD[lev]->VyMP,               parH[lev]->memSizerealkMP,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->VzMP,          parD[lev]->VzMP,               parH[lev]->memSizerealkMP,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->RhoMP,         parD[lev]->RhoMP,              parH[lev]->memSizerealkMP,  cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeMeasurePointsIndex(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->kMP));
	checkCudaErrors( cudaFreeHost(parH[lev]->VxMP));
	checkCudaErrors( cudaFreeHost(parH[lev]->VyMP));
	checkCudaErrors( cudaFreeHost(parH[lev]->VzMP));
	checkCudaErrors( cudaFreeHost(parH[lev]->RhoMP));
}
void Parameter::cudaAllocFsForCheckPointAndRestart(int lev)
{
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->d0SP.f[0] ),           (unsigned long long)getD3Qxx()*(unsigned long long)parH[lev]->mem_size_real_SP));
}
void Parameter::cudaCopyFsForRestart(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->d0SP.f[0],  parH[lev]->d0SP.f[0],     (unsigned long long)getD3Qxx()*(unsigned long long)parH[lev]->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyFsForCheckPoint(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->d0SP.f[0],  parD[lev]->d0SP.f[0],     (unsigned long long)getD3Qxx()*(unsigned long long)parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeFsForCheckPointAndRestart(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->d0SP.f[0]));
}
//DragLift
void Parameter::cudaAllocDragLift(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(double)*numofelem;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPreX), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPreY), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPreZ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPostX), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPostY), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->DragPostZ), mem_size  ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPreX), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPreY), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPreZ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPostX), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPostY), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->DragPostZ), mem_size  ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 6. * (double)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyDragLift(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(double)*numofelem;

	checkCudaErrors( cudaMemcpy(parH[lev]->DragPreX, parD[lev]->DragPreX, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->DragPreY, parD[lev]->DragPreY, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->DragPreZ, parD[lev]->DragPreZ, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->DragPostX, parD[lev]->DragPostX, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->DragPostY, parD[lev]->DragPostY, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->DragPostZ, parD[lev]->DragPostZ, mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeDragLift(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPreX));
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPreY));
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPreZ));
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPostX));
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPostY));
	checkCudaErrors( cudaFreeHost(parH[lev]->DragPostZ));
}
//2ndMoments
void Parameter::cudaAlloc2ndMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kxyFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kyzFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kxzFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kxxMyyFromfcNEQ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->kxxMzzFromfcNEQ), mem_size  ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kxyFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kyzFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kxzFromfcNEQ   ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kxxMyyFromfcNEQ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->kxxMzzFromfcNEQ), mem_size  ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 5. * (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopy2ndMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	checkCudaErrors( cudaMemcpy(parH[lev]->kxyFromfcNEQ   , parD[lev]->kxyFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->kyzFromfcNEQ   , parD[lev]->kyzFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->kxzFromfcNEQ   , parD[lev]->kxzFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->kxxMyyFromfcNEQ, parD[lev]->kxxMyyFromfcNEQ, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->kxxMzzFromfcNEQ, parD[lev]->kxxMzzFromfcNEQ, mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFree2ndMoments(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->kxyFromfcNEQ   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->kyzFromfcNEQ   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->kxzFromfcNEQ   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->kxxMyyFromfcNEQ));
	checkCudaErrors( cudaFreeHost(parH[lev]->kxxMzzFromfcNEQ));
}
//3rdMoments
void Parameter::cudaAlloc3rdMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbbb ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMabc ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbac ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbca ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcba ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMacb ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcab ), mem_size ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbbb ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMabc ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbac ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbca ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcba ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMacb ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcab ), mem_size ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 7. * (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopy3rdMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbbb, parD[lev]->CUMbbb, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMabc, parD[lev]->CUMabc, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbac, parD[lev]->CUMbac, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbca, parD[lev]->CUMbca, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcba, parD[lev]->CUMcba, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMacb, parD[lev]->CUMacb, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcab, parD[lev]->CUMcab, mem_size, cudaMemcpyDeviceToHost));
}																														   
void Parameter::cudaFree3rdMoments(int lev)																				   
{
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbbb ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMabc ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbac ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbca ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcba ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMacb ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcab ));
}
//higher order moments
void Parameter::cudaAllocHigherMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcbb ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbcb ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbbc ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcca ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcac ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMacc ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMbcc ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMcbc ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMccb ), mem_size ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->CUMccc ), mem_size ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcbb ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbcb ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbbc ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcca ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcac ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMacc ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMbcc ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMcbc ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMccb ), mem_size ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->CUMccc ), mem_size ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 7. * (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyHigherMoments(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcbb, parD[lev]->CUMcbb, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbcb, parD[lev]->CUMbcb, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbbc, parD[lev]->CUMbbc, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcca, parD[lev]->CUMcca, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcac, parD[lev]->CUMcac, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMacc, parD[lev]->CUMacc, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMbcc, parD[lev]->CUMbcc, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMcbc, parD[lev]->CUMcbc, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMccb, parD[lev]->CUMccb, mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->CUMccc, parD[lev]->CUMccc, mem_size, cudaMemcpyDeviceToHost));
}																														   
void Parameter::cudaFreeHigherMoments(int lev)																				   
{
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcbb ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbcb ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbbc ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcca ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcac ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMacc ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMbcc ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMcbc ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMccb ));
	checkCudaErrors( cudaFreeHost(parH[lev]->CUMccc ));
}
//Velcities to fit the Forcing
void Parameter::cudaAllocForceVelo(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VxForce   ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VyForce   ), mem_size  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->VzForce   ), mem_size  ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VxForce   ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VyForce   ), mem_size  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->VzForce   ), mem_size  ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyForceVelo(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;

	checkCudaErrors( cudaMemcpy(parH[lev]->VxForce   , parD[lev]->VxForce   , mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->VyForce   , parD[lev]->VyForce   , mem_size, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->VzForce   , parD[lev]->VzForce   , mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeForceVelo(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->VxForce   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->VyForce   ));
	checkCudaErrors( cudaFreeHost(parH[lev]->VzForce   ));
}
//Forcing
void Parameter::cudaAllocForcing()
{
	unsigned int mem_size = sizeof(real) * 3;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(this->forcingH ), mem_size));
	//Device
	checkCudaErrors( cudaMalloc(    (void**) &(this->forcingD ), mem_size));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyForcingToDevice()
{
	unsigned int mem_size = sizeof(real) * 3;
	checkCudaErrors( cudaMemcpy(this->forcingD, this->forcingH , mem_size, cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyForcingToHost()
{
	unsigned int mem_size = sizeof(real) * 3;
	checkCudaErrors( cudaMemcpy(this->forcingH, this->forcingD , mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeForcing()
{
	checkCudaErrors( cudaFreeHost(this->forcingH));
}
//cp Top
void Parameter::cudaAllocCpTop(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpTop;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpTop;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpPressTop), mem_size_double  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpTopIndex), mem_size_int     ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpPressTop), mem_size_double      ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpTopIndex), mem_size_int         ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_double + (double)mem_size_int;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyCpTopInit(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpTop;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpTop;

	checkCudaErrors( cudaMemcpy(parD[lev]->cpPressTop, parH[lev]->cpPressTop, mem_size_double, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->cpTopIndex, parH[lev]->cpTopIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyCpTop(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpTop;
	//unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpTop;

	checkCudaErrors( cudaMemcpy(parH[lev]->cpPressTop, parD[lev]->cpPressTop, mem_size_double, cudaMemcpyDeviceToHost));
	//checkCudaErrors( cudaMemcpy(parH[lev]->cpTopIndex, parD[lev]->cpTopIndex, mem_size_int,    cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeCpTop(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->cpPressTop));
	checkCudaErrors( cudaFreeHost(parH[lev]->cpTopIndex));
}
//cp Bottom
void Parameter::cudaAllocCpBottom(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpBottom;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpPressBottom), mem_size_double  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpBottomIndex), mem_size_int     ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpPressBottom), mem_size_double      ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpBottomIndex), mem_size_int         ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_double + (double)mem_size_int;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyCpBottomInit(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpBottom;

	checkCudaErrors( cudaMemcpy(parD[lev]->cpPressBottom, parH[lev]->cpPressBottom, mem_size_double, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->cpBottomIndex, parH[lev]->cpBottomIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyCpBottom(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom;
	//unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpBottom;

	checkCudaErrors( cudaMemcpy(parH[lev]->cpPressBottom, parD[lev]->cpPressBottom, mem_size_double, cudaMemcpyDeviceToHost));
	//checkCudaErrors( cudaMemcpy(parH[lev]->cpBottomIndex, parD[lev]->cpBottomIndex, mem_size_int,    cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeCpBottom(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->cpPressBottom));
	checkCudaErrors( cudaFreeHost(parH[lev]->cpBottomIndex));
}
//cp Bottom 2
void Parameter::cudaAllocCpBottom2(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom2;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpBottom2;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpPressBottom2), mem_size_double  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->cpBottom2Index), mem_size_int     ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpPressBottom2), mem_size_double      ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->cpBottom2Index), mem_size_int         ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_double + (double)mem_size_int;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyCpBottom2Init(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom2;
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsCpBottom2;

	checkCudaErrors( cudaMemcpy(parD[lev]->cpPressBottom2, parH[lev]->cpPressBottom2, mem_size_double, cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->cpBottom2Index, parH[lev]->cpBottom2Index, mem_size_int,    cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyCpBottom2(int lev)
{
	unsigned int mem_size_double = sizeof(double)       * parH[lev]->numberOfPointsCpBottom2;

	checkCudaErrors( cudaMemcpy(parH[lev]->cpPressBottom2, parD[lev]->cpPressBottom2, mem_size_double, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeCpBottom2(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->cpPressBottom2));
	checkCudaErrors( cudaFreeHost(parH[lev]->cpBottom2Index));
}
//////////////////////////////////////////////////////////////////////////
//particles
void Parameter::cudaAllocParticles(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordXlocal),        parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordYlocal),        parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordZlocal),        parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordXabsolut),      parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordYabsolut),      parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.coordZabsolut),      parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.veloX),              parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.veloY),              parH[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.veloZ),              parH[lev]->plp.memSizerealAll  ));
	//checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.randomLocationInit), parH[lev]->plp.memSizereal     ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.ID),                 parH[lev]->plp.memSizeID          ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.cellBaseID),         parH[lev]->plp.memSizeID          ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.timestep),           parH[lev]->plp.memSizeTimestep    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.stuck),              parH[lev]->plp.memSizeBool        ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->plp.hot),                parH[lev]->plp.memSizeBoolBC      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordXlocal),            parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordYlocal),            parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordZlocal),            parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordXabsolut),          parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordYabsolut),          parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.coordZabsolut),          parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.veloX),                  parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.veloY),                  parD[lev]->plp.memSizerealAll  ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.veloZ),                  parD[lev]->plp.memSizerealAll  ));
	//checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.randomLocationInit),     parD[lev]->plp.memSizereal     ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.ID),                     parD[lev]->plp.memSizeID          ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.cellBaseID),             parD[lev]->plp.memSizeID          ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.timestep),               parD[lev]->plp.memSizeTimestep    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.stuck),                  parD[lev]->plp.memSizeBool        ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->plp.hot),                    parD[lev]->plp.memSizeBoolBC      ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)parD[lev]->plp.memSizerealAll * (double)9.0 + (double)parD[lev]->plp.memSizeID * (double)2.0 + (double)parD[lev]->plp.memSizeTimestep 
		+ (double)parD[lev]->plp.memSizeBool + (double)parD[lev]->plp.memSizeBoolBC;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyParticles(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordXlocal,        parD[lev]->plp.coordXlocal,        parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordYlocal,        parD[lev]->plp.coordYlocal,        parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordZlocal,        parD[lev]->plp.coordZlocal,        parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordXabsolut,      parD[lev]->plp.coordXabsolut,      parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordYabsolut,      parD[lev]->plp.coordYabsolut,      parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.coordZabsolut,      parD[lev]->plp.coordZabsolut,      parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.veloX,              parD[lev]->plp.veloX,              parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.veloY,              parD[lev]->plp.veloY,              parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.veloZ,              parD[lev]->plp.veloZ,              parH[lev]->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
	//checkCudaErrors( cudaMemcpy(parH[lev]->plp.randomLocationInit, parD[lev]->plp.randomLocationInit, parH[lev]->plp.memSizereal,     cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.ID,                 parD[lev]->plp.ID,                 parH[lev]->plp.memSizeID,          cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.cellBaseID,         parD[lev]->plp.cellBaseID,         parH[lev]->plp.memSizeID,          cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parH[lev]->plp.timestep,           parD[lev]->plp.timestep,           parH[lev]->plp.memSizeTimestep,    cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeParticles(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordXlocal)       );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordYlocal)       );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordZlocal)       );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordXabsolut)     );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordYabsolut)     );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.coordZabsolut)     );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.veloX)             );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.veloY)             );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.veloZ)             );
	//checkCudaErrors( cudaFreeHost(parH[lev]->plp.randomLocationInit));
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.ID)                );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.cellBaseID)        );
	checkCudaErrors( cudaFreeHost(parH[lev]->plp.timestep)          );
}
//random values
void Parameter::cudaAllocRandomValues()
{
	//Device
	checkCudaErrors( cudaMalloc((void**)&(this->devState), (sizeof(curandState)*parD[getFine()]->plp.numberOfParticles) ));
}
//////////////////////////////////////////////////////////////////////////
//porous media
void Parameter::cudaAllocPorousMedia(PorousMedia* pm, int lev)
{
	unsigned int mem_size_IDsPM = sizeof(unsigned int)*pm->getSizePM();
	unsigned int *tmpIDHost, *tmpIDDevice;
	//std::cout << "cudaMallocHost" << endl;
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(tmpIDHost), mem_size_IDsPM));

	//std::cout << "cudaMalloc" << endl;
	//Device
	checkCudaErrors(cudaMalloc((void**) &(tmpIDDevice), mem_size_IDsPM));

	//std::cout << "set Host and Device arrays PM" << endl;
	//////////////////////////////////////////////////////////////////////////
	pm->setHostNodeIDsPM(tmpIDHost);
	pm->setDeviceNodeIDsPM(tmpIDDevice);
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_IDsPM;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPorousMedia(PorousMedia* pm, int lev)
{
	unsigned int mem_size_IDsPM = sizeof(unsigned int)*pm->getSizePM();
	unsigned int *tmpIDHost   = pm->getHostNodeIDsPM();
	unsigned int *tmpIDDevice = pm->getDeviceNodeIDsPM();
	//////////////////////////////////////////////////////////////////////////
	checkCudaErrors(cudaMemcpy(tmpIDDevice, tmpIDHost, mem_size_IDsPM, cudaMemcpyHostToDevice));
	//////////////////////////////////////////////////////////////////////////
	pm->setDeviceNodeIDsPM(tmpIDDevice);
}
void Parameter::cudaFreePorousMedia(PorousMedia* pm, int lev)
{
	checkCudaErrors(cudaFreeHost(pm->getHostNodeIDsPM()));
}
//////////////////////////////////////////////////////////////////////////
//advection diffusion
void Parameter::cudaAllocConc(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->Conc), parH[lev]->mem_size_real_SP));	
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->Conc), parD[lev]->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyConcDH(int lev)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->Conc, parD[lev]->Conc,  parH[lev]->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void Parameter::cudaCopyConcHD(int lev)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->Conc, parH[lev]->Conc, parH[lev]->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeConc(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->Conc));
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocTempFs(int lev)
{
	//Device
	if (getDiffMod() == 7)
	{
		checkCudaErrors( cudaMalloc((void**) &(parD[lev]->d7.f[0]), getDiffMod()*parH[lev]->mem_size_real_SP));
	} 
	else if (getDiffMod() == 27)
	{
		checkCudaErrors( cudaMalloc((void**) &(parD[lev]->d27.f[0]), getDiffMod()*parH[lev]->mem_size_real_SP));
	}	
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)(getDiffMod()*parH[lev]->mem_size_real_SP);
	setMemsizeGPU(tmp, false);
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocTempPressBC(int lev)
{
	unsigned int mem_size_TempPress_k = sizeof(int)*parH[lev]->TempPress.kTemp;
	unsigned int mem_size_TempPress_q = sizeof(real)*parH[lev]->TempPress.kTemp;

	// Host Memory
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempPress.temp, mem_size_TempPress_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempPress.velo, mem_size_TempPress_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempPress.k,    mem_size_TempPress_k ));

	// Device Memory
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempPress.temp, mem_size_TempPress_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempPress.velo, mem_size_TempPress_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempPress.k,    mem_size_TempPress_k));

}
void Parameter::cudaCopyTempPressBCHD(int lev)
{
	unsigned int mem_size_TempPress_k = sizeof(int)*parH[lev]->TempPress.kTemp;
	unsigned int mem_size_TempPress_q = sizeof(real)*parH[lev]->TempPress.kTemp;

	checkCudaErrors( cudaMemcpy(parD[lev]->TempPress.temp, parH[lev]->TempPress.temp, mem_size_TempPress_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->TempPress.velo, parH[lev]->TempPress.velo, mem_size_TempPress_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->TempPress.k,    parH[lev]->TempPress.k,    mem_size_TempPress_k,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeTempPressBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->TempPress.temp));
	checkCudaErrors( cudaFreeHost(parH[lev]->TempPress.velo));
	checkCudaErrors( cudaFreeHost(parH[lev]->TempPress.k   ));
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocTempVeloBC(int lev)
{
	unsigned int mem_size_TempVel_k = sizeof(int)*parH[lev]->TempVel.kTemp;
	unsigned int mem_size_TempVel_q = sizeof(real)*parH[lev]->TempVel.kTemp;

	printf("mem_size_TempVel_k = %d,  mem_size_TempVel_q = %d \n", mem_size_TempVel_k, mem_size_TempVel_q);
	// Host Memory
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempVel.temp,      mem_size_TempVel_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempVel.tempPulse, mem_size_TempVel_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempVel.velo,      mem_size_TempVel_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->TempVel.k,         mem_size_TempVel_k ));

	// Device Memory
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempVel.temp,      mem_size_TempVel_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempVel.tempPulse, mem_size_TempVel_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempVel.velo,      mem_size_TempVel_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->TempVel.k,         mem_size_TempVel_k));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)(mem_size_TempVel_q * 3.0 + mem_size_TempVel_k);
	setMemsizeGPU(tmp, false);

}
void Parameter::cudaCopyTempVeloBCHD(int lev)
{
	unsigned int mem_size_TempVel_k = sizeof(int)*parH[lev]->TempVel.kTemp;
	unsigned int mem_size_TempVel_q = sizeof(real)*parH[lev]->TempVel.kTemp;

	printf("mem_size_TempVel_k = %d,  mem_size_TempVel_q = %d \n", mem_size_TempVel_k, mem_size_TempVel_q);
	checkCudaErrors( cudaMemcpy(parD[lev]->TempVel.temp,      parH[lev]->TempVel.temp,      mem_size_TempVel_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->TempVel.tempPulse, parH[lev]->TempVel.tempPulse, mem_size_TempVel_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->TempVel.velo,      parH[lev]->TempVel.velo,      mem_size_TempVel_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->TempVel.k,         parH[lev]->TempVel.k,         mem_size_TempVel_k,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeTempVeloBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->TempVel.temp     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->TempVel.tempPulse));
	checkCudaErrors( cudaFreeHost(parH[lev]->TempVel.velo     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->TempVel.k        ));
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocTempNoSlipBC(int lev)
{
	unsigned int mem_size_Temp_k = sizeof(int)*parH[lev]->Temp.kTemp;
	unsigned int mem_size_Temp_q = sizeof(real)*parH[lev]->Temp.kTemp;

	// Host Memory
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->Temp.temp, mem_size_Temp_q ));
	checkCudaErrors( cudaMallocHost((void**) &parH[lev]->Temp.k,    mem_size_Temp_k ));

	// Device Memory
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->Temp.temp, mem_size_Temp_q));
	checkCudaErrors( cudaMalloc((void**) &parD[lev]->Temp.k,    mem_size_Temp_k));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)(mem_size_Temp_q + mem_size_Temp_k);
	setMemsizeGPU(tmp, false);

}
void Parameter::cudaCopyTempNoSlipBCHD(int lev)
{
	unsigned int mem_size_Temp_k = sizeof(int)*parH[lev]->Temp.kTemp;
	unsigned int mem_size_Temp_q = sizeof(real)*parH[lev]->Temp.kTemp;

	checkCudaErrors( cudaMemcpy(parD[lev]->Temp.temp, parH[lev]->Temp.temp, mem_size_Temp_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parD[lev]->Temp.k,    parH[lev]->Temp.k,    mem_size_Temp_k,  cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeTempNoSlipBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->Temp.temp));
	checkCudaErrors( cudaFreeHost(parH[lev]->Temp.k   ));
}
//PlaneConc
void Parameter::cudaAllocPlaneConcIn(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->ConcPlaneIn), mem_size  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->ConcPlaneIn), mem_size  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPlaneConcIn(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	checkCudaErrors( cudaMemcpy(parH[lev]->ConcPlaneIn,   parD[lev]->ConcPlaneIn,   mem_size, cudaMemcpyDeviceToHost));
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocPlaneConcOut1(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->ConcPlaneOut1), mem_size  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->ConcPlaneOut1), mem_size  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPlaneConcOut1(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	checkCudaErrors( cudaMemcpy(parH[lev]->ConcPlaneOut1, parD[lev]->ConcPlaneOut1, mem_size, cudaMemcpyDeviceToHost));
}
//////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocPlaneConcOut2(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->ConcPlaneOut2), mem_size  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->ConcPlaneOut2), mem_size  ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyPlaneConcOut2(int lev, int numofelem)
{
	unsigned int mem_size = sizeof(real)*numofelem;
	checkCudaErrors( cudaMemcpy(parH[lev]->ConcPlaneOut2, parD[lev]->ConcPlaneOut2, mem_size, cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreePlaneConc(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->ConcPlaneIn));
	checkCudaErrors( cudaFreeHost(parH[lev]->ConcPlaneOut1));
	checkCudaErrors( cudaFreeHost(parH[lev]->ConcPlaneOut2));
}
//////////////////////////////////////////////////////////////////////////
//concentration file
void Parameter::cudaAllocConcFile(int lev)
{
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsConc;

	printf("numberOfPointsConc = %d \n", mem_size_int);
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->concIndex), mem_size_int     ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->concIndex), mem_size_int         ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_int;
	setMemsizeGPU(tmp, false);
}
void Parameter::cudaCopyConcFile(int lev)
{
	unsigned int mem_size_int    = sizeof(unsigned int) * parH[lev]->numberOfPointsConc;

	checkCudaErrors( cudaMemcpy(parD[lev]->concIndex, parH[lev]->concIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void Parameter::cudaFreeConcFile(int lev)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->concIndex));
}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//Process Neighbors
//1D domain decomposition
void Parameter::cudaAllocProcessNeighbor(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighbor[processNeighbor].index ),                  parH[lev]->sendProcessNeighbor[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighbor[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->sendProcessNeighbor[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighbor[processNeighbor].index ),                  parH[lev]->recvProcessNeighbor[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighbor[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->recvProcessNeighbor[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighbor[processNeighbor].index ),                      parD[lev]->sendProcessNeighbor[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighbor[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->sendProcessNeighbor[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighbor[processNeighbor].index ),                      parD[lev]->recvProcessNeighbor[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighbor[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->recvProcessNeighbor[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighbor[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->sendProcessNeighbor[processNeighbor].memsizeFs + 
				 (double)parH[lev]->recvProcessNeighbor[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->recvProcessNeighbor[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighbor[processNeighbor].index, 
								parH[lev]->sendProcessNeighbor[processNeighbor].index, 
								parH[lev]->sendProcessNeighbor[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighbor[processNeighbor].index, 
								parH[lev]->recvProcessNeighbor[processNeighbor].index, 
								parH[lev]->recvProcessNeighbor[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighbor[processNeighbor].f[0], 
								parH[lev]->recvProcessNeighbor[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->recvProcessNeighbor[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighbor[processNeighbor].f[0], 
								parD[lev]->sendProcessNeighbor[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->sendProcessNeighbor[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighbor(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighbor[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighbor[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighbor[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighbor[processNeighbor].f[0]     ));
}
////////////////////////////////////////////////////////////////////////////////////
//  3D domain decomposition
//  X  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborX[processNeighbor].index ),                  parH[lev]->sendProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborX[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->sendProcessNeighborX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborX[processNeighbor].index ),                  parH[lev]->recvProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborX[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->recvProcessNeighborX[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborX[processNeighbor].index ),                      parD[lev]->sendProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborX[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->sendProcessNeighborX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborX[processNeighbor].index ),                      parD[lev]->recvProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborX[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->recvProcessNeighborX[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborX[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->sendProcessNeighborX[processNeighbor].memsizeFs + 
				 (double)parH[lev]->recvProcessNeighborX[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->recvProcessNeighborX[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborX[processNeighbor].index, 
								parH[lev]->sendProcessNeighborX[processNeighbor].index, 
								parH[lev]->sendProcessNeighborX[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborX[processNeighbor].index, 
								parH[lev]->recvProcessNeighborX[processNeighbor].index, 
								parH[lev]->recvProcessNeighborX[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborX[processNeighbor].f[0], 
								parH[lev]->recvProcessNeighborX[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->recvProcessNeighborX[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborX[processNeighbor].f[0], 
								parD[lev]->sendProcessNeighborX[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->sendProcessNeighborX[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborX[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborX[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborX[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborX[processNeighbor].f[0]     ));
}
//  Y  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborY[processNeighbor].index ),                  parH[lev]->sendProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborY[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->sendProcessNeighborY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborY[processNeighbor].index ),                  parH[lev]->recvProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborY[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->recvProcessNeighborY[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborY[processNeighbor].index ),                      parD[lev]->sendProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborY[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->sendProcessNeighborY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborY[processNeighbor].index ),                      parD[lev]->recvProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborY[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->recvProcessNeighborY[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborY[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->sendProcessNeighborY[processNeighbor].memsizeFs + 
				 (double)parH[lev]->recvProcessNeighborY[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->recvProcessNeighborY[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborY[processNeighbor].index, 
								parH[lev]->sendProcessNeighborY[processNeighbor].index, 
								parH[lev]->sendProcessNeighborY[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborY[processNeighbor].index, 
								parH[lev]->recvProcessNeighborY[processNeighbor].index, 
								parH[lev]->recvProcessNeighborY[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborY[processNeighbor].f[0], 
								parH[lev]->recvProcessNeighborY[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->recvProcessNeighborY[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborY[processNeighbor].f[0], 
								parD[lev]->sendProcessNeighborY[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->sendProcessNeighborY[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborY[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborY[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborY[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborY[processNeighbor].f[0]     ));
}
//  Z  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborZ[processNeighbor].index ),                  parH[lev]->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborZ[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborZ[processNeighbor].index ),                  parH[lev]->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborZ[processNeighbor].f[0]  ),     getD3Qxx() * parH[lev]->recvProcessNeighborZ[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborZ[processNeighbor].index ),                      parD[lev]->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborZ[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborZ[processNeighbor].index ),                      parD[lev]->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborZ[processNeighbor].f[0]  ),         getD3Qxx() * parD[lev]->recvProcessNeighborZ[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborZ[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->sendProcessNeighborZ[processNeighbor].memsizeFs + 
				 (double)parH[lev]->recvProcessNeighborZ[processNeighbor].memsizeIndex + (double)getD3Qxx()*(double)parH[lev]->recvProcessNeighborZ[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborZ[processNeighbor].index, 
								parH[lev]->sendProcessNeighborZ[processNeighbor].index, 
								parH[lev]->sendProcessNeighborZ[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborZ[processNeighbor].index, 
								parH[lev]->recvProcessNeighborZ[processNeighbor].index, 
								parH[lev]->recvProcessNeighborZ[processNeighbor].memsizeIndex, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborZ[processNeighbor].f[0], 
								parH[lev]->recvProcessNeighborZ[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->recvProcessNeighborZ[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborZ[processNeighbor].f[0], 
								parD[lev]->sendProcessNeighborZ[processNeighbor].f[0], 
								getD3Qxx() * parD[lev]->sendProcessNeighborZ[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborZ[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborZ[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborZ[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborZ[processNeighbor].f[0]     ));
}
////////////////////////////////////////////////////////////////////////////////////
//  3D domain decomposition convection diffusion
//  X  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborADX(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADX[processNeighbor].index ),                  parH[lev]->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADX[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADX[processNeighbor].index ),                  parH[lev]->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADX[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADX[processNeighbor].index ),                      parD[lev]->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADX[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADX[processNeighbor].index ),                      parD[lev]->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADX[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborADX[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->sendProcessNeighborADX[processNeighbor].memsizeFs + 
		         (double)parH[lev]->recvProcessNeighborADX[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->recvProcessNeighborADX[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborADXIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborADX[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADX[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADX[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADX[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADX[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADX[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADXFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADX[processNeighbor].f[0], 
		                        parH[lev]->recvProcessNeighborADX[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->recvProcessNeighborADX[processNeighbor].memsizeFs, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADXFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborADX[processNeighbor].f[0], 
		                        parD[lev]->sendProcessNeighborADX[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->sendProcessNeighborADX[processNeighbor].memsizeFs, 
		                        cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborADX(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADX[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADX[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADX[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADX[processNeighbor].f[0]     ));
}
//  Y  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborADY(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADY[processNeighbor].index ),                  parH[lev]->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADY[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADY[processNeighbor].index ),                  parH[lev]->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADY[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADY[processNeighbor].index ),                      parD[lev]->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADY[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADY[processNeighbor].index ),                      parD[lev]->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADY[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborADY[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->sendProcessNeighborADY[processNeighbor].memsizeFs + 
		         (double)parH[lev]->recvProcessNeighborADY[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->recvProcessNeighborADY[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborADYIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborADY[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADY[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADY[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADY[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADY[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADY[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADYFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADY[processNeighbor].f[0], 
		                        parH[lev]->recvProcessNeighborADY[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->recvProcessNeighborADY[processNeighbor].memsizeFs, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADYFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborADY[processNeighbor].f[0], 
		                        parD[lev]->sendProcessNeighborADY[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->sendProcessNeighborADY[processNeighbor].memsizeFs, 
		                        cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborADY(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADY[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADY[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADY[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADY[processNeighbor].f[0]     ));
}
//  Z  /////////////////////////////////////////////////////////////////////////////
void Parameter::cudaAllocProcessNeighborADZ(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADZ[processNeighbor].index ),                  parH[lev]->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->sendProcessNeighborADZ[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADZ[processNeighbor].index ),                  parH[lev]->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parH[lev]->recvProcessNeighborADZ[processNeighbor].f[0]  ),   getDiffMod() * parH[lev]->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADZ[processNeighbor].index ),                      parD[lev]->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->sendProcessNeighborADZ[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADZ[processNeighbor].index ),                      parD[lev]->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parD[lev]->recvProcessNeighborADZ[processNeighbor].f[0]  ),       getDiffMod() * parD[lev]->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp = (double)parH[lev]->sendProcessNeighborADZ[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->sendProcessNeighborADZ[processNeighbor].memsizeFs + 
		         (double)parH[lev]->recvProcessNeighborADZ[processNeighbor].memsizeIndex + (double)getDiffMod()*(double)parH[lev]->recvProcessNeighborADZ[processNeighbor].memsizeFs;
	setMemsizeGPU(tmp, false);
	//printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void Parameter::cudaCopyProcessNeighborADZIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors( cudaMemcpy(parD[lev]->sendProcessNeighborADZ[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADZ[processNeighbor].index, 
		                        parH[lev]->sendProcessNeighborADZ[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADZ[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADZ[processNeighbor].index, 
		                        parH[lev]->recvProcessNeighborADZ[processNeighbor].memsizeIndex, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADZFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parD[lev]->recvProcessNeighborADZ[processNeighbor].f[0], 
		                        parH[lev]->recvProcessNeighborADZ[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->recvProcessNeighborADZ[processNeighbor].memsizeFs, 
		                        cudaMemcpyHostToDevice));
}
void Parameter::cudaCopyProcessNeighborADZFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parH[lev]->sendProcessNeighborADZ[processNeighbor].f[0], 
		                        parD[lev]->sendProcessNeighborADZ[processNeighbor].f[0], 
		         getDiffMod() * parD[lev]->sendProcessNeighborADZ[processNeighbor].memsizeFs, 
		                        cudaMemcpyDeviceToHost));
}
void Parameter::cudaFreeProcessNeighborADZ(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADZ[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parH[lev]->sendProcessNeighborADZ[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADZ[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parH[lev]->recvProcessNeighborADZ[processNeighbor].f[0]     ));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//set-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::setForcing(real forcingX, real forcingY, real forcingZ)
{
	this->hostForcing[0] = forcingX;
	this->hostForcing[1] = forcingY;
	this->hostForcing[2] = forcingZ;
}
void Parameter::setPhi(real inPhi)
{
	Phi = inPhi;
}
void Parameter::setAngularVelocity(real inAngVel)
{
	angularVelocity = inAngVel;
}
void Parameter::setStepEnsight(unsigned int step)
{
	this->stepEnsight = step;
}
void Parameter::setOutputCount(unsigned int outputCount)
{
	this->outputCount = outputCount;
}
void Parameter::setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK)
{
	this->limitOfNodesForVTK = limitOfNodesForVTK;
}
void Parameter::setStartTurn(unsigned int inStartTurn)
{
	startTurn = inStartTurn;
}
void Parameter::setDiffOn(bool isDiff)
{
	diffOn = isDiff;
}
void Parameter::setDiffMod(int DiffMod)
{
	diffMod = DiffMod;
}
void Parameter::setD3Qxx(int d3qxx)
{
	this->D3Qxx = d3qxx;
}
void Parameter::setMaxLevel(int maxlevel)
{
	this->maxlevel = maxlevel-1;
}
void Parameter::setParticleBasicLevel(int pbl)
{
	this->particleBasicLevel = pbl;
}
void Parameter::setParticleInitLevel(int pil)
{
	this->particleInitLevel = pil;
}
void Parameter::setNumberOfParticles(int nop)
{
	this->numberOfParticles = nop;
}
void Parameter::setCalcParticles(bool calcParticles)
{
	this->calcParticles = calcParticles;
}
void Parameter::setStartXHotWall(real startXHotWall)
{
	this->startXHotWall = startXHotWall;
}
void Parameter::setEndXHotWall(real endXHotWall)
{
	this->endXHotWall = endXHotWall;
}
void Parameter::setTEnd(unsigned int tend)
{
	ic.tend = tend;
}
void Parameter::setTOut(unsigned int tout)
{
	ic.tout = tout;
}
void Parameter::setTStartOut(unsigned int tStartOut)
{
	ic.tStartOut = tStartOut;
}
void Parameter::setCalcMedian(bool calcMedian)
{
	ic.calcMedian = calcMedian;
}
void Parameter::setTimeCalcMedStart(int CalcMedStart)
{
	ic.tCalcMedStart = CalcMedStart;
}
void Parameter::setTimeCalcMedEnd(int CalcMedEnd)
{
	ic.tCalcMedEnd = CalcMedEnd;
}
void Parameter::setOutputPath(std::string oPath)
{
	ic.oPath = oPath;
}
void Parameter::setOutputPrefix(std::string oPrefix)
{
	//std::string test = fname;
	ic.oPrefix = oPrefix;
}
void Parameter::setFName(std::string fname)
{
	//std::string test = fname;
	ic.fname = fname;
}
void Parameter::setPrintFiles(bool printfiles)
{
	ic.printFiles = printfiles;
}
void Parameter::setReadGeo(bool readGeo)
{
	ic.readGeo = readGeo;
}
void Parameter::setDiffusivity(real Diffusivity)
{
	ic.Diffusivity = Diffusivity;
}
void Parameter::setTemperatureInit(real Temp)
{
	ic.Temp = Temp;
}
void Parameter::setTemperatureBC(real TempBC)
{
	ic.TempBC = TempBC;
}
void Parameter::setViscosity(real Viscosity)
{
	ic.vis = Viscosity;
}
void Parameter::setVelocity(real Velocity)
{
	ic.u0 = Velocity;
}
void Parameter::setViscosityRatio(real ViscosityRatio)
{
	ic.vis_ratio = ViscosityRatio;
}
void Parameter::setVelocityRatio(real VelocityRatio)
{
	ic.u0_ratio = VelocityRatio;
}
void Parameter::setDensityRatio(real DensityRatio)
{
	ic.delta_rho = DensityRatio;
}
void Parameter::setPressRatio(real PressRatio)
{
	ic.delta_press = PressRatio;
}
void Parameter::setRealX(real RealX)
{
	ic.RealX = RealX;
}
void Parameter::setRealY(real RealY)
{
	ic.RealY = RealY;
}
void Parameter::setPressInID(unsigned int PressInID)
{
	ic.PressInID = PressInID;
}
void Parameter::setPressOutID(unsigned int PressOutID)
{
	ic.PressOutID = PressOutID;
}
void Parameter::setPressInZ(unsigned int PressInZ)
{
	ic.PressInZ = PressInZ;
}
void Parameter::setPressOutZ(unsigned int PressOutZ)
{
	ic.PressOutZ = PressOutZ;
}
void Parameter::setMaxDev(int maxdev)
{
	ic.maxdev = maxdev;
}
void Parameter::setMyID(int myid)
{
	ic.myid = myid;
}
void Parameter::setNumprocs(int numprocs)
{
	ic.numprocs = numprocs;
}
void Parameter::setDevices(std::vector<int> devices)
{
	ic.devices = devices;
}
void Parameter::setGeometryFileC(std::string GeometryFileC)
{
	ic.geometryFileC = GeometryFileC;
}
void Parameter::setGeometryFileM(std::string GeometryFileM)
{
	ic.geometryFileM = GeometryFileM;
}
void Parameter::setGeometryFileF(std::string GeometryFileF)
{
	ic.geometryFileF = GeometryFileF;
}
void Parameter::setNeedInterface(std::vector<bool> NeedInterface)
{
	ic.NeedInterface = NeedInterface;
}
void Parameter::setRe(real Re)
{
	ic.Re = Re;
}
void Parameter::setFactorPressBC(real factorPressBC)
{
	ic.factorPressBC = factorPressBC;
}
void Parameter::setIsGeo(bool isGeo)
{
	ic.isGeo = isGeo;
}
void Parameter::setIsGeoNormal(bool isGeoNormal)
{
	ic.isGeoNormal = isGeoNormal;
}
void Parameter::setIsInflowNormal(bool isInflowNormal)
{
	ic.isInflowNormal = isInflowNormal;
}
void Parameter::setIsOutflowNormal(bool isOutflowNormal)
{
	ic.isOutflowNormal = isOutflowNormal;
}
void Parameter::setIsProp(bool isProp)
{
	ic.isProp = isProp;
}
void Parameter::setIsCp(bool isCp)
{
	ic.isCp = isCp;
}
void Parameter::setConcFile(bool concFile)
{
	ic.isConc = concFile;
}
void Parameter::setUseMeasurePoints(bool useMeasurePoints)
{
	ic.isMeasurePoints = useMeasurePoints;
}
void Parameter::setUseWale(bool useWale)
{
	ic.isWale = useWale;
}
void Parameter::setSimulatePorousMedia(bool simulatePorousMedia)
{
	ic.simulatePorousMedia = simulatePorousMedia;
}
void Parameter::setGridX(std::vector<int> GridX)
{
	ic.GridX = GridX;
}
void Parameter::setGridY(std::vector<int> GridY)
{
	ic.GridY = GridY;
}
void Parameter::setGridZ(std::vector<int> GridZ)
{
	ic.GridZ = GridZ;
}
void Parameter::setDistX(std::vector<int> DistX)
{
	ic.DistX = DistX;
}
void Parameter::setDistY(std::vector<int> DistY)
{
	ic.DistY = DistY;
}
void Parameter::setDistZ(std::vector<int> DistZ)
{
	ic.DistZ = DistZ;
}
void Parameter::setScaleLBMtoSI(std::vector<real> scaleLBMtoSI)
{
	ic.scaleLBMtoSI = scaleLBMtoSI;
}
void Parameter::setTranslateLBMtoSI(std::vector<real> translateLBMtoSI)
{
	ic.translateLBMtoSI = translateLBMtoSI;
}
void Parameter::setMinCoordX(std::vector<real> MinCoordX)
{
	ic.minCoordX = MinCoordX;
}
void Parameter::setMinCoordY(std::vector<real> MinCoordY)
{
	ic.minCoordY = MinCoordY;
}
void Parameter::setMinCoordZ(std::vector<real> MinCoordZ)
{
	ic.minCoordZ = MinCoordZ;
}
void Parameter::setMaxCoordX(std::vector<real> MaxCoordX)
{
	ic.maxCoordX = MaxCoordX;
}
void Parameter::setMaxCoordY(std::vector<real> MaxCoordY)
{
	ic.maxCoordY = MaxCoordY;
}
void Parameter::setMaxCoordZ(std::vector<real> MaxCoordZ)
{
	ic.maxCoordZ = MaxCoordZ;
}
void Parameter::setTempH(TempforBoundaryConditions* TempH)
{
	this->TempH = TempH;
}
void Parameter::setTempD(TempforBoundaryConditions* TempD)
{
	this->TempD = TempD;
}
void Parameter::setTempVelH(TempVelforBoundaryConditions* TempVelH)
{
	this->TempVelH = TempVelH;
}
void Parameter::setTempVelD(TempVelforBoundaryConditions* TempVelD)
{
	this->TempVelD = TempVelD;
}
void Parameter::setTempPressH(TempPressforBoundaryConditions* TempPressH)
{
	this->TempPressH = TempPressH;
}
void Parameter::setTempPressD(TempPressforBoundaryConditions* TempPressD)
{
	this->TempPressD = TempPressD;
}
//void Parameter::setkInflowQ(unsigned int kInflowQ)
//{
//   this->kInflowQ = kInflowQ;
//}
//void Parameter::setkOutflowQ(unsigned int kOutflowQ)
//{
//   this->kOutflowQ = kOutflowQ;
//}
//void Parameter::setQinflowH(QforBoundaryConditions* QinflowH)
//{
//   this->QinflowH = QinflowH;
//}
//void Parameter::setQinflowD(QforBoundaryConditions* QinflowD)
//{
//   this->QinflowD = QinflowD;
//}
//void Parameter::setQoutflowH(QforBoundaryConditions* QoutflowH)
//{
//   this->QoutflowH = QoutflowH;
//}
//void Parameter::setQoutflowD(QforBoundaryConditions* QoutflowD)
//{
//   this->QoutflowD = QoutflowD;
//}
void Parameter::setkFull(std::string kFull)
{
	ic.kFull = kFull;
}
void Parameter::setgeoFull(std::string geoFull)
{
	ic.geoFull = geoFull;
}
void Parameter::setgeoVec(std::string geoVec)
{
	ic.geoVec = geoVec;
}
void Parameter::setcoordX(std::string coordX)
{
	ic.coordX = coordX;
}
void Parameter::setcoordY(std::string coordY)
{
	ic.coordY = coordY;
}
void Parameter::setcoordZ(std::string coordZ)
{
	ic.coordZ = coordZ;
}
void Parameter::setneighborX(std::string neighborX)
{
	ic.neighborX = neighborX;
}
void Parameter::setneighborY(std::string neighborY)
{
	ic.neighborY = neighborY;
}
void Parameter::setneighborZ(std::string neighborZ)
{
	ic.neighborZ = neighborZ;
}
void Parameter::setneighborWSB(std::string neighborWSB)
{
	ic.neighborWSB = neighborWSB;
}
void Parameter::setscaleCFC(std::string scaleCFC)
{
	ic.scaleCFC = scaleCFC;
}
void Parameter::setscaleCFF(std::string scaleCFF)
{
	ic.scaleCFF = scaleCFF;
}
void Parameter::setscaleFCC(std::string scaleFCC)
{
	ic.scaleFCC = scaleFCC;
}
void Parameter::setscaleFCF(std::string scaleFCF)
{
	ic.scaleFCF = scaleFCF;
}
void Parameter::setscaleOffsetCF(std::string scaleOffsetCF)
{
	ic.scaleOffsetCF = scaleOffsetCF;
}
void Parameter::setscaleOffsetFC(std::string scaleOffsetFC)
{
	ic.scaleOffsetFC = scaleOffsetFC;
}
void Parameter::setgeomBoundaryBcQs(std::string geomBoundaryBcQs)
{
	ic.geomBoundaryBcQs = geomBoundaryBcQs;
}
void Parameter::setgeomBoundaryBcValues(std::string geomBoundaryBcValues)
{
	ic.geomBoundaryBcValues = geomBoundaryBcValues;
}
void Parameter::setnoSlipBcPos(std::string noSlipBcPos)
{
	ic.noSlipBcPos = noSlipBcPos;
}
void Parameter::setnoSlipBcQs(std::string noSlipBcQs)
{
	ic.noSlipBcQs = noSlipBcQs;
}
void Parameter::setnoSlipBcValue(std::string noSlipBcValue)
{
	ic.noSlipBcValue = noSlipBcValue;
}
void Parameter::setnoSlipBcValues(std::string noSlipBcValues)
{
	ic.noSlipBcValues = noSlipBcValues;
}
void Parameter::setslipBcPos(std::string slipBcPos)
{
	ic.slipBcPos = slipBcPos;
}
void Parameter::setslipBcQs(std::string slipBcQs)
{
	ic.slipBcQs = slipBcQs;
}
void Parameter::setslipBcValue(std::string slipBcValue)
{
	ic.slipBcValue = slipBcValue;
}
void Parameter::setpressBcPos(std::string pressBcPos)
{
	ic.pressBcPos = pressBcPos;
}
void Parameter::setpressBcQs(std::string pressBcQs)
{
	ic.pressBcQs = pressBcQs;
}
void Parameter::setpressBcValue(std::string pressBcValue)
{
	ic.pressBcValue = pressBcValue;
}
void Parameter::setpressBcValues(std::string pressBcValues)
{
	ic.pressBcValues = pressBcValues;
}
void Parameter::setvelBcQs(std::string velBcQs)
{
	ic.velBcQs = velBcQs;
}
void Parameter::setvelBcValues(std::string velBcValues)
{
	ic.velBcValues = velBcValues;
}
void Parameter::setinletBcQs(std::string inletBcQs)
{
	ic.inletBcQs = inletBcQs;
}
void Parameter::setinletBcValues(std::string inletBcValues)
{
	ic.inletBcValues = inletBcValues;
}
void Parameter::setoutletBcQs(std::string outletBcQs)
{
	ic.outletBcQs = outletBcQs;
}
void Parameter::setoutletBcValues(std::string outletBcValues)
{
	ic.outletBcValues = outletBcValues;
}
void Parameter::settopBcQs(std::string topBcQs)
{
	ic.topBcQs = topBcQs;
}
void Parameter::settopBcValues(std::string topBcValues)
{
	ic.topBcValues = topBcValues;
}
void Parameter::setbottomBcQs(std::string bottomBcQs)
{
	ic.bottomBcQs = bottomBcQs;
}
void Parameter::setbottomBcValues(std::string bottomBcValues)
{
	ic.bottomBcValues = bottomBcValues;
}
void Parameter::setfrontBcQs(std::string frontBcQs)
{
	ic.frontBcQs = frontBcQs;
}
void Parameter::setfrontBcValues(std::string frontBcValues)
{
	ic.frontBcValues = frontBcValues;
}
void Parameter::setbackBcQs(std::string backBcQs)
{
	ic.backBcQs = backBcQs;
}
void Parameter::setbackBcValues(std::string backBcValues)
{
	ic.backBcValues = backBcValues;
}
void Parameter::setwallBcQs(std::string wallBcQs)
{
	ic.wallBcQs = wallBcQs;
}
void Parameter::setwallBcValues(std::string wallBcValues)
{
	ic.wallBcValues = wallBcValues;
}
void Parameter::setperiodicBcQs(std::string periodicBcQs)
{
	ic.periodicBcQs = periodicBcQs;
}
void Parameter::setperiodicBcValues(std::string periodicBcValues)
{
	ic.periodicBcValues = periodicBcValues;
}
void Parameter::setpropellerQs(std::string propellerQs)
{
	ic.propellerQs = propellerQs;
}
void Parameter::setpropellerValues(std::string propellerValues)
{
	ic.propellerValues = propellerValues;
}
void Parameter::setpropellerCylinder(std::string propellerCylinder)
{
	ic.propellerCylinder = propellerCylinder;
}
void Parameter::setmeasurePoints(std::string measurePoints)
{
	ic.measurePoints = measurePoints;
}
void Parameter::setnumberNodes(std::string numberNodes)
{
	ic.numberNodes = numberNodes;
}
void Parameter::setLBMvsSI(std::string LBMvsSI)
{
	ic.LBMvsSI = LBMvsSI;
}
void Parameter::setcpTop(std::string cpTop)
{
	ic.cpTop = cpTop;
}
void Parameter::setcpBottom(std::string cpBottom)
{
	ic.cpBottom = cpBottom;
}
void Parameter::setcpBottom2(std::string cpBottom2)
{
	ic.cpBottom2 = cpBottom2;
}
void Parameter::setConcentration(std::string concFile)
{
	ic.concentration = concFile;
}
void Parameter::setclockCycleForMP(real clockCycleForMP)
{
	ic.clockCycleForMP = clockCycleForMP;
}
void Parameter::setTimeDoCheckPoint(unsigned int tDoCheckPoint)
{
	ic.tDoCheckPoint = tDoCheckPoint;
}
void Parameter::setTimeDoRestart(unsigned int tDoRestart)
{
	ic.tDoRestart = tDoRestart;
}
void Parameter::setDoCheckPoint(bool doCheckPoint)
{
	ic.doCheckPoint = doCheckPoint;
}
void Parameter::setDoRestart(bool doRestart)
{
	ic.doRestart = doRestart;
}
void Parameter::settimestepForMP(unsigned int timestepForMP)
{
	ic.timeStepForMP = timestepForMP;
}
void Parameter::setObj(std::string str, bool isObj)
{
	if (str == "geo")
	{
		this->setIsGeo(isObj);
	}
	else if (str == "prop")
	{
		this->setIsProp(isObj);
	}
	else if (str == "cp")
	{
		this->setIsCp(isObj);
	}
	else if (str == "geoNormal")
	{
		this->setIsGeoNormal(isObj);
	}
	else if (str == "inflowNormal")
	{
		this->setIsInflowNormal(isObj);
	}
	else if (str == "outflowNormal")
	{
		this->setIsOutflowNormal(isObj);
	}
}
void Parameter::setGeometryValues(bool GeometryValues)
{
	ic.GeometryValues = GeometryValues;
}
void Parameter::setCalc2ndOrderMoments(bool is2ndOrderMoments)
{
	ic.is2ndOrderMoments = is2ndOrderMoments;
}
void Parameter::setCalc3rdOrderMoments(bool is3rdOrderMoments)
{
	ic.is3rdOrderMoments = is3rdOrderMoments;
}
void Parameter::setCalcHighOrderMoments(bool isHighOrderMoments)
{
	ic.isHighOrderMoments = isHighOrderMoments;
}
void Parameter::setMemsizeGPU(double admem, bool reset)
{
	if (reset == true)
	{
		this->memsizeGPU = 0.;
	} 
	else
	{
		this->memsizeGPU += admem;
	}
}
//1D domain decomposition
void Parameter::setPossNeighborFiles(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSend = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecv = possNeighborFiles;
	}
}
void Parameter::setNumberOfProcessNeighbors(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
	}
}
void Parameter::setIsNeighbor(bool isNeigbor)
{
	this->isNeigbor = isNeigbor;
}
//3D domain decomposition
void Parameter::setPossNeighborFilesX(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendX = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvX = possNeighborFiles;
	}
}
void Parameter::setPossNeighborFilesY(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendY = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvY = possNeighborFiles;
	}
}
void Parameter::setPossNeighborFilesZ(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendZ = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvZ = possNeighborFiles;
	}
}
void Parameter::setNumberOfProcessNeighborsX(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
			parD[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	}
}
void Parameter::setNumberOfProcessNeighborsY(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
			parD[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	}
}
void Parameter::setNumberOfProcessNeighborsZ(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->recvProcessNeighborADZ.resize(numberOfProcessNeighbors);
			parD[level]->recvProcessNeighborADZ.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	}
}
void Parameter::setIsNeighborX(bool isNeigbor)
{
	this->isNeigborX = isNeigbor;
}
void Parameter::setIsNeighborY(bool isNeigbor)
{
	this->isNeigborY = isNeigbor;
}
void Parameter::setIsNeighborZ(bool isNeigbor)
{
	this->isNeigborZ = isNeigbor;
}
void Parameter::setgeomBoundaryNormalX(std::string geomNormalX)
{
	ic.geomNormalX = geomNormalX;
}
void Parameter::setgeomBoundaryNormalY(std::string geomNormalY)
{
	ic.geomNormalY = geomNormalY;
}
void Parameter::setgeomBoundaryNormalZ(std::string geomNormalZ)
{
	ic.geomNormalZ = geomNormalZ;
}
void Parameter::setInflowBoundaryNormalX(std::string inflowNormalX)
{
	ic.inflowNormalX = inflowNormalX;
}
void Parameter::setInflowBoundaryNormalY(std::string inflowNormalY)
{
	ic.inflowNormalY = inflowNormalY;
}
void Parameter::setInflowBoundaryNormalZ(std::string inflowNormalZ)
{
	ic.inflowNormalZ = inflowNormalZ;
}
void Parameter::setOutflowBoundaryNormalX(std::string outflowNormalX)
{
	ic.outflowNormalX = outflowNormalX;
}
void Parameter::setOutflowBoundaryNormalY(std::string outflowNormalY)
{
	ic.outflowNormalY = outflowNormalY;
}
void Parameter::setOutflowBoundaryNormalZ(std::string outflowNormalZ)
{
	ic.outflowNormalZ = outflowNormalZ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double* Parameter::getForcesDouble()
{
	return this->hostForcing;
}
real* Parameter::getForcesHost()
{
	return this->forcingH;
}
real* Parameter::getForcesDev()
{
	return this->forcingD;
}
real Parameter::getPhi()
{
	return Phi;
}
real Parameter::getAngularVelocity()
{
	return angularVelocity;
}
real Parameter::getStartXHotWall()
{
	return this->startXHotWall;
}
real Parameter::getEndXHotWall()
{
	return this->endXHotWall;
}
unsigned int Parameter::getStepEnsight()
{
	return this->stepEnsight;
}
unsigned int Parameter::getOutputCount()
{
	return this->outputCount;
}
unsigned int Parameter::getlimitOfNodesForVTK()
{
	return this->limitOfNodesForVTK;
}
unsigned int Parameter::getStartTurn()
{
	return startTurn;
}
ParameterStruct* Parameter::getParD(int level)
{
	return parD[level];
}
ParameterStruct* Parameter::getParH(int level)
{
	return parH[level];
}
unsigned int Parameter::getSizeMat(int level)
{
	return parH[level]->size_Mat;
}
unsigned int Parameter::getMemSizereal(int level)
{
	return parH[level]->mem_size_real;
}
unsigned int Parameter::getMemSizeInt(int level)
{
	return parH[level]->mem_size_int;
}
unsigned int Parameter::getMemSizeBool(int level)
{
	return parH[level]->mem_size_bool;
}
unsigned int Parameter::getMemSizerealYZ(int level)
{
	return parH[level]->mem_size_real_yz;
}
int Parameter::getFine()
{
	return fine;
}
int Parameter::getCoarse()
{
	return coarse;
}
int Parameter::getParticleBasicLevel()
{
	return this->particleBasicLevel;
}
int Parameter::getParticleInitLevel()
{
	return this->particleInitLevel;
}
int Parameter::getNumberOfParticles()
{
	return this->numberOfParticles;
}
bool Parameter::getEvenOrOdd(int level)
{
	return parH[level]->evenOrOdd;
}
bool Parameter::getDiffOn()
{
	return diffOn;
}
int Parameter::getDiffMod()
{
	return diffMod;
}
int Parameter::getFactorNZ()
{
	return factor_gridNZ;
}
int Parameter::getD3Qxx()
{
	return this->D3Qxx;
}
int Parameter::getMaxLevel()
{
	return this->maxlevel;
}
unsigned int Parameter::getTStart()
{
	if (getDoRestart())
	{
		return getTimeDoRestart() + 1;
	} 
	else
	{
		return 1;
	}
}
unsigned int Parameter::getTInit()
{
	if (getDoRestart())
	{
		return getTimeDoRestart();
	} 
	else
	{
		return 0;
	}
}
unsigned int Parameter::getTEnd()
{
	return ic.tend;
}
unsigned int Parameter::getTOut()
{
	return ic.tout;
}
unsigned int Parameter::getTStartOut()
{
	return ic.tStartOut;
}
bool Parameter::getCalcMedian()
{
	return ic.calcMedian;
}
bool Parameter::getCalcParticle()
{
	return this->calcParticles;
}
int Parameter::getTimeCalcMedStart()
{
	return ic.tCalcMedStart;
}
int Parameter::getTimeCalcMedEnd()
{
	return ic.tCalcMedEnd;
}
std::string Parameter::getOutputPath()
{
	return ic.oPath;
}
std::string Parameter::getOutputPrefix()
{
	return ic.oPrefix;
}
std::string Parameter::getFName()
{
	return ic.fname;
}
bool Parameter::getPrintFiles()
{
	return ic.printFiles;
}
bool Parameter::getReadGeo()
{
	return ic.readGeo;
}
real Parameter::getDiffusivity()
{
	return ic.Diffusivity;
}
real Parameter::getTemperatureInit()
{
	return ic.Temp;
}
real Parameter::getTemperatureBC()
{
	return ic.TempBC;
}
real Parameter::getViscosity()
{
	return ic.vis;
}
real Parameter::getVelocity()
{
	return ic.u0;
}
real Parameter::getViscosityRatio()
{
	return ic.vis_ratio;
}
real Parameter::getVelocityRatio()
{
	return ic.u0_ratio;
}
real Parameter::getDensityRatio()
{
	return ic.delta_rho;
}
real Parameter::getPressRatio()
{
	return ic.delta_press;
}
real Parameter::getRealX()
{
	return ic.RealX;
}
real Parameter::getRealY()
{
	return ic.RealY;
}
unsigned int Parameter::getPressInID()
{
	return ic.PressInID;
}
unsigned int Parameter::getPressOutID()
{
	return ic.PressOutID;
}
unsigned int Parameter::getPressInZ()
{
	return ic.PressInZ;
}
unsigned int Parameter::getPressOutZ()
{
	return ic.PressOutZ;
}
int Parameter::getMaxDev()
{
	return ic.maxdev;
}
int Parameter::getMyID()
{
	return ic.myid;
}
int Parameter::getNumprocs()
{
	return ic.numprocs;
}
std::vector<int> Parameter::getDevices()
{
	return ic.devices;
}
std::string Parameter::getGeometryFileC()
{
	return ic.geometryFileC;
}
std::string Parameter::getGeometryFileM()
{
	return ic.geometryFileM;
}
std::string Parameter::getGeometryFileF()
{
	return ic.geometryFileF;
}
std::vector<bool> Parameter::getNeedInterface()
{
	return ic.NeedInterface;
}
real Parameter::getRe()
{
	return ic.Re;
}
real Parameter::getFactorPressBC()
{
	return ic.factorPressBC;
}
std::vector<int> Parameter::getGridX()
{
	return ic.GridX;
}
std::vector<int> Parameter::getGridY()
{
	return ic.GridY;
}
std::vector<int> Parameter::getGridZ()
{
	return ic.GridZ;
}
std::vector<int> Parameter::getDistX()
{
	return ic.DistX;
}
std::vector<int> Parameter::getDistY()
{
	return ic.DistY;
}
std::vector<int> Parameter::getDistZ()
{
	return ic.DistZ;
}
std::vector<real> Parameter::getScaleLBMtoSI()
{
	return ic.scaleLBMtoSI;
}
std::vector<real> Parameter::getTranslateLBMtoSI()
{
	return ic.translateLBMtoSI;
}
std::vector<real> Parameter::getMinCoordX()
{
	return ic.minCoordX;
}
std::vector<real> Parameter::getMinCoordY()
{
	return ic.minCoordY;
}
std::vector<real> Parameter::getMinCoordZ()
{
	return ic.minCoordZ;
}
std::vector<real> Parameter::getMaxCoordX()
{
	return ic.maxCoordX;
}
std::vector<real> Parameter::getMaxCoordY()
{
	return ic.maxCoordY;
}
std::vector<real> Parameter::getMaxCoordZ()
{
	return ic.maxCoordZ;
}
TempforBoundaryConditions* Parameter::getTempH()
{
	return this->TempH;
}
TempforBoundaryConditions* Parameter::getTempD()
{
	return this->TempD;
}
TempVelforBoundaryConditions* Parameter::getTempVelH()
{
	return this->TempVelH;
}
TempVelforBoundaryConditions* Parameter::getTempVelD()
{
	return this->TempVelD;
}
TempPressforBoundaryConditions* Parameter::getTempPressH()
{
	return this->TempPressH;
}
TempPressforBoundaryConditions* Parameter::getTempPressD()
{
	return this->TempPressD;
}
//unsigned int Parameter::getkInflowQ()
//{
//   return this->kInflowQ;
//}
//unsigned int Parameter::getkOutflowQ()
//{
//   return this->kOutflowQ;
//}
//QforBoundaryConditions* Parameter::getQinflowH()
//{
//   return this->QinflowH;
//}
//QforBoundaryConditions* Parameter::getQinflowD()
//{
//   return this->QinflowD;
//}
//QforBoundaryConditions* Parameter::getQoutflowH()
//{
//   return this->QoutflowH;
//}
//QforBoundaryConditions* Parameter::getQoutflowD()
//{
//   return this->QoutflowD;
//}
std::string Parameter::getkFull()
{
	return ic.kFull;
}
std::string Parameter::getgeoFull()
{
	return ic.geoFull;
}
std::string Parameter::getgeoVec()
{
	return ic.geoVec;
}
std::string Parameter::getcoordX()
{
	return ic.coordX;
}
std::string Parameter::getcoordY()
{
	return ic.coordY;
}
std::string Parameter::getcoordZ()
{
	return ic.coordZ;
}
std::string Parameter::getneighborX()
{
	return ic.neighborX;
}
std::string Parameter::getneighborY()
{
	return ic.neighborY;
}
std::string Parameter::getneighborZ()
{
	return ic.neighborZ;
}
std::string Parameter::getneighborWSB()
{
	return ic.neighborWSB;
}
std::string Parameter::getscaleCFC()
{
	return ic.scaleCFC;
}
std::string Parameter::getscaleCFF()
{
	return ic.scaleCFF;
}
std::string Parameter::getscaleFCC()
{
	return ic.scaleFCC;
}
std::string Parameter::getscaleFCF()
{
	return ic.scaleFCF;
}
std::string Parameter::getscaleOffsetCF()
{
	return ic.scaleOffsetCF;
}
std::string Parameter::getscaleOffsetFC()
{
	return ic.scaleOffsetFC;
}
std::string Parameter::getgeomBoundaryBcQs()
{
	return ic.geomBoundaryBcQs;
}
std::string Parameter::getgeomBoundaryBcValues()
{
	return ic.geomBoundaryBcValues;
}
std::string Parameter::getnoSlipBcPos()
{
	return ic.noSlipBcPos;
}
std::string Parameter::getnoSlipBcQs()
{
	return ic.noSlipBcQs;
}
std::string Parameter::getnoSlipBcValue()
{
	return ic.noSlipBcValue;
}
std::string Parameter::getnoSlipBcValues()
{
	return ic.noSlipBcValues;
}
std::string Parameter::getslipBcPos()
{
	return ic.slipBcPos;
}
std::string Parameter::getslipBcQs()
{
	return ic.slipBcQs;
}
std::string Parameter::getslipBcValue()
{
	return ic.slipBcValue;
}
std::string Parameter::getpressBcPos()
{
	return ic.pressBcPos;
}
std::string Parameter::getpressBcQs()
{
	return ic.pressBcQs;
}
std::string Parameter::getpressBcValue()
{
	return ic.pressBcValue;
}
std::string Parameter::getpressBcValues()
{
	return ic.pressBcValues;
}
std::string Parameter::getvelBcQs()
{
	return ic.velBcQs;
}
std::string Parameter::getvelBcValues()
{
	return ic.velBcValues;
}
std::string Parameter::getinletBcQs()
{
	return ic.inletBcQs;
}
std::string Parameter::getinletBcValues()
{
	return ic.inletBcValues;
}
std::string Parameter::getoutletBcQs()
{
	return ic.outletBcQs;
}
std::string Parameter::getoutletBcValues()
{
	return ic.outletBcValues;
}
std::string Parameter::gettopBcQs()
{
	return ic.topBcQs;
}
std::string Parameter::gettopBcValues()
{
	return ic.topBcValues;
}
std::string Parameter::getbottomBcQs()
{
	return ic.bottomBcQs;
}
std::string Parameter::getbottomBcValues()
{
	return ic.bottomBcValues;
}
std::string Parameter::getfrontBcQs()
{
	return ic.frontBcQs;
}
std::string Parameter::getfrontBcValues()
{
	return ic.frontBcValues;
}
std::string Parameter::getbackBcQs()
{
	return ic.backBcQs;
}
std::string Parameter::getbackBcValues()
{
	return ic.backBcValues;
}
std::string Parameter::getwallBcQs()
{
	return ic.wallBcQs;
}
std::string Parameter::getwallBcValues()
{
	return ic.wallBcValues;
}
std::string Parameter::getperiodicBcQs()
{
	return ic.periodicBcQs;
}
std::string Parameter::getperiodicBcValues()
{
	return ic.periodicBcValues;
}
std::string Parameter::getpropellerQs()
{
	return ic.propellerQs;
}
std::string Parameter::getpropellerValues()
{
	return ic.propellerValues;
}
std::string Parameter::getpropellerCylinder()
{
	return ic.propellerCylinder;
}
std::string Parameter::getmeasurePoints()
{
	return ic.measurePoints;
}
std::string Parameter::getLBMvsSI()
{
	return ic.LBMvsSI;
}
std::string Parameter::getnumberNodes()
{
	return ic.numberNodes;
}
std::string Parameter::getcpTop()
{
	return ic.cpTop;
}
std::string Parameter::getcpBottom()
{
	return ic.cpBottom;
}
std::string Parameter::getcpBottom2()
{
	return ic.cpBottom2;
}
std::string Parameter::getConcentration()
{
	return ic.concentration;
}
real Parameter::getclockCycleForMP()
{
	return ic.clockCycleForMP;
}
unsigned int Parameter::getTimeDoCheckPoint()
{
	return ic.tDoCheckPoint;
}
unsigned int Parameter::getTimeDoRestart()
{
	return ic.tDoRestart;
}
bool Parameter::getDoCheckPoint()
{
	return ic.doCheckPoint;
}
bool Parameter::getDoRestart()
{
	return ic.doRestart;
}
bool Parameter::getIsGeo()
{
	return ic.isGeo;
}
bool Parameter::getIsGeoNormal()
{
	return ic.isGeoNormal;
}
bool Parameter::getIsInflowNormal()
{
	return ic.isInflowNormal;
}
bool Parameter::getIsOutflowNormal()
{
	return ic.isOutflowNormal;
}
bool Parameter::getIsCp()
{
	return ic.isCp;
}
bool Parameter::getConcFile()
{
	return ic.isConc;
}
bool Parameter::getUseMeasurePoints()
{
	return ic.isMeasurePoints;
}
bool Parameter::getUseWale()
{
	return ic.isWale;
}
bool Parameter::getSimulatePorousMedia()
{
	return ic.simulatePorousMedia;
}
bool Parameter::getIsGeometryValues()
{
	return ic.GeometryValues;
}
bool Parameter::getCalc2ndOrderMoments()
{
	return ic.is2ndOrderMoments;
}
bool Parameter::getCalc3rdOrderMoments()
{
	return ic.is3rdOrderMoments;
}
bool Parameter::getCalcHighOrderMoments()
{
	return ic.isHighOrderMoments;
}
bool Parameter::getIsProp()
{
	return ic.isProp;
}
bool Parameter::overWritingRestart(unsigned int t)
{
	if (t == getTimeDoRestart())
	{
		return true;
	} 
	else
	{
		return false;
	}
}
unsigned int Parameter::getTimestepForMP()
{
	return ic.timeStepForMP;
}
double Parameter::getMemsizeGPU()
{
	return this->memsizeGPU;
}
//1D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFiles(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSend;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecv;
	}
}
unsigned int Parameter::getNumberOfProcessNeighbors(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighbor.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighbor.size();
	}
}
bool Parameter::getIsNeighbor()
{
	return this->isNeigbor;
}
//3D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFilesX(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendX;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvX;
	}
}
std::vector<std::string> Parameter::getPossNeighborFilesY(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendY;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvY;
	}
}
std::vector<std::string> Parameter::getPossNeighborFilesZ(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendZ;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvZ;
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsX(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborX.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborX.size();
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsY(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborY.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborY.size();
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsZ(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborZ.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborZ.size();
	}
}
bool Parameter::getIsNeighborX()
{
	return this->isNeigborX;
}
bool Parameter::getIsNeighborY()
{
	return this->isNeigborY;
}
bool Parameter::getIsNeighborZ()
{
	return this->isNeigborZ;
}
std::string Parameter::getgeomBoundaryNormalX()
{
	return ic.geomNormalX;
}
std::string Parameter::getgeomBoundaryNormalY()
{
	return ic.geomNormalY;
}
std::string Parameter::getgeomBoundaryNormalZ()
{
	return ic.geomNormalZ;
}
std::string Parameter::getInflowBoundaryNormalX()
{
	return ic.inflowNormalX;
}
std::string Parameter::getInflowBoundaryNormalY()
{
	return ic.inflowNormalY;
}
std::string Parameter::getInflowBoundaryNormalZ()
{
	return ic.inflowNormalZ;
}
std::string Parameter::getOutflowBoundaryNormalX()
{
	return ic.outflowNormalX;
}
std::string Parameter::getOutflowBoundaryNormalY()
{
	return ic.outflowNormalY;
}
std::string Parameter::getOutflowBoundaryNormalZ()
{
	return ic.outflowNormalZ;
}
curandState* Parameter::getRandomState()
{
	return this->devState;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::initInterfaceParameter(int level)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//host
	parH[level]->K_CF               = 0;
	parH[level]->K_FC               = 0;
	if (parH[level]->need_interface[INTERFACE_E]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2);
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2);
	} 
	if (parH[level]->need_interface[INTERFACE_W]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_N]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_S]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_T]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_B]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2);
	}
	//parH[level]->K_CF               = (( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2)*2)+
	//                                  (( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2)*2)+
	//                                  (( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2)*2);
	//parH[level]->K_FC               = (((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2)*2)+
	//                                  (((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2)*2)+
	//                                  (((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2)*2);
	parH[level]->mem_size_kCF       = sizeof(unsigned int)*parH[level]->K_CF;
	parH[level]->mem_size_kFC       = sizeof(unsigned int)*parH[level]->K_FC;
	parH[level]->mem_size_kCF_off   = sizeof(real)*parH[level]->K_CF;
	parH[level]->mem_size_kFC_off   = sizeof(real)*parH[level]->K_FC;
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//device
	parD[level]->K_CF               = parH[level]->K_CF;
	parD[level]->K_FC               = parH[level]->K_FC;
	parD[level]->mem_size_kCF       = parH[level]->mem_size_kCF;
	parD[level]->mem_size_kFC       = parH[level]->mem_size_kFC;
	parD[level]->mem_size_kCF_off   = parH[level]->mem_size_kCF_off;
	parD[level]->mem_size_kFC_off   = parH[level]->mem_size_kFC_off;
	///////////////////////////////////////////////////////////////////////////////////////////////////
}
real Parameter::TrafoXtoWorld(int CoordX, int level)
{
	return (parH[level]->mTtoWx*CoordX+parH[level]->cTtoWx);
}
real Parameter::TrafoYtoWorld(int CoordY, int level)
{
	return (parH[level]->mTtoWy*CoordY+parH[level]->cTtoWy);
}
real Parameter::TrafoZtoWorld(int CoordZ, int level)
{
	return (parH[level]->mTtoWz*CoordZ+parH[level]->cTtoWz);
}
real Parameter::TrafoXtoMGsWorld(int CoordX, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->XdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordX ) * parH[level]->dx);
	return temp;
}
real Parameter::TrafoYtoMGsWorld(int CoordY, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->YdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordY ) * parH[level]->dx);
	return temp;
}
real Parameter::TrafoZtoMGsWorld(int CoordZ, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->ZdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordZ) * parH[level]->dx);
	return temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
