#include "Interface.h"
#include <algorithm>
#include <math.h>
#define TEST false



Interface::Interface(bool binaer)
{
	this->binaer=binaer;
	system = new string[6];
	system[0] = "inlet";
	system[1] = "outlet";
	system[2] = "front";
	system[3] = "back";
	system[4] = "top";
	system[5] = "bottom";



}

Interface::~Interface(void)
{
}

bool Interface::getBinaer(){
	return binaer;
}

void Interface::allocArrays_CoordNeighborGeo(Parameter* para) {
	cout << "-----Config Arrays Coord, Neighbor, Geo------" << endl;

	CoordNeighborGeoV *coordX = new CoordNeighborGeoV(para->getcoordX(), binaer, true);
	CoordNeighborGeoV *coordY = new CoordNeighborGeoV(para->getcoordY(), binaer, true);
	CoordNeighborGeoV *coordZ = new CoordNeighborGeoV(para->getcoordZ(), binaer, true);
	neighX = new CoordNeighborGeoV(para->getneighborX(), binaer, false);
	neighY = new CoordNeighborGeoV(para->getneighborY(), binaer, false);
	neighZ = new CoordNeighborGeoV(para->getneighborZ(), binaer, false);
	//Particles or Wale
	if (para->getCalcParticle() || para->getUseWale())
	{
		neighWSB = new CoordNeighborGeoV(para->getneighborWSB(), binaer, false);
	}

	CoordNeighborGeoV *geoV = new CoordNeighborGeoV(para->getgeoVec(), binaer, false);


	int level = coordX->getLevel();
	cout << "Anzahl Level: " << level + 1 << endl;
	int AnzahlKnotenGes = 0;
	cout << "Anzahl Knoten: " << endl;
	//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------
	for (int i = 0; i <= level; i++) {
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		int temp = coordX->getSize(i) + 1;
		AnzahlKnotenGes += temp;
		cout << "Level " << i << " = " << temp << " Knoten" << endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int blocks = (temp / para->getParH(i)->numberofthreads) + 1;

		//cout << "Anzahl Blocks vorher: " << temp / para->getParH(i)->numberofthreads << endl;
		//cout << "Anzahl Blocks nachher: " << blocks << endl;

		para->getParH(i)->size_Array_SP = blocks * para->getParH(i)->numberofthreads;
		para->getParD(i)->size_Array_SP = para->getParH(i)->size_Array_SP;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->getParH(i)->size_Mat_SP = temp;
		para->getParD(i)->size_Mat_SP = temp;
		para->getParH(i)->mem_size_doubflo_SP = sizeof(doubflo)* para->getParH(i)->size_Array_SP;
		para->getParH(i)->mem_size_int_SP = sizeof(unsigned int)* para->getParH(i)->size_Array_SP;
		para->getParD(i)->mem_size_doubflo_SP = sizeof(doubflo)* para->getParD(i)->size_Array_SP;
		para->getParD(i)->mem_size_int_SP = sizeof(unsigned int)* para->getParD(i)->size_Array_SP;
		//para->getParH(i)->mem_size_doubflo_SP = sizeof(doubflo)* para->getParH(i)->size_Mat_SP;
		//para->getParH(i)->mem_size_int_SP = sizeof(unsigned int)* para->getParH(i)->size_Mat_SP;
		//para->getParD(i)->mem_size_doubflo_SP = sizeof(doubflo)* para->getParD(i)->size_Mat_SP;
		//para->getParD(i)->mem_size_int_SP = sizeof(unsigned int)* para->getParD(i)->size_Mat_SP;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//cout <<"Test 1 " <<endl;
		//cout << "size_Mat_SP Host in Level " << i << " = " << para->getParH(i)->size_Mat_SP << endl;
		//cout << "size_Mat_SP Device in Level " << i << " = " << para->getParD(i)->size_Mat_SP << endl;
		//cout << "mem_size_doubflo_SP Host in Level " << i << " = " << para->getParH(i)->mem_size_doubflo_SP << endl;
		//cout << "mem_size_doubflo_SP Device in Level " << i << " = " << para->getParD(i)->mem_size_doubflo_SP << endl;
		//cout << "mem_size_int_SP Host in Level " << i << " = " << para->getParH(i)->mem_size_int_SP << endl;
		//cout << "mem_size_int_SP Device in Level " << i << " = " << para->getParD(i)->mem_size_int_SP << endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		para->cudaAllocCoord(i);
		para->cudaAllocSP(i);
		if (para->getCalcMedian()) 	                       
			para->cudaAllocMedianSP(i);
		if (para->getCalcParticle() || para->getUseWale()) 
			para->cudaAllocNeighborWSB(i);
		if (para->getUseWale())                            
			para->cudaAllocTurbulentViscosity(i);
		//cout <<"Test 2 " <<endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		coordX->initArrayCoord(para->getParH(i)->coordX_SP, i);
		coordY->initArrayCoord(para->getParH(i)->coordY_SP, i);
		coordZ->initArrayCoord(para->getParH(i)->coordZ_SP, i);
		//cout <<"Test 3 " <<endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		neighX->initArray(para->getParH(i)->neighborX_SP, i);
		neighY->initArray(para->getParH(i)->neighborY_SP, i);
		neighZ->initArray(para->getParH(i)->neighborZ_SP, i);
		//Particles or Wale
		if (para->getCalcParticle() || para->getUseWale())
		{
			neighWSB->initArray(para->getParH(i)->neighborWSB_SP, i);
		}
		//cout <<"Test 4 " <<endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		geoV->initArray(para->getParH(i)->geoSP, i);
		//cout <<"Test 5 " <<endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////test
		//cout << " Laenge des Vektors: " << coordX1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << coordY1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << coordZ1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << neighX1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << neighY1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << neighZ1->getSize(i) << endl;
		//cout << " Laenge des Vektors: " << geoV1->getSize(i) << endl;
		////////////////////////////////////////////////////////////////////////////////
		//Test Coarse Nodes in Fine Grid
		//if (i == 0)
		//{
		//	for (int t1 = 0; t1 < para->getParH(0)->size_Mat_SP; t1++)
		//	{
		//		if (((para->getParH(0)->coordX_SP[t1] > 15.0) && (para->getParH(0)->coordX_SP[t1] < 25.0)) &&
		//			((para->getParH(0)->coordY_SP[t1] > 15.0) && (para->getParH(0)->coordY_SP[t1] < 25.0)) &&
		//			((para->getParH(0)->coordZ_SP[t1] > 15.0) && (para->getParH(0)->coordZ_SP[t1] < 25.0)))
		//		{
		//			cout << " CoordTest " << para->getParH(0)->coordX_SP[t1] << " " << para->getParH(0)->coordY_SP[t1] << " " << para->getParH(0)->coordZ_SP[t1] << " " << t1 << endl;
		//		}
		//		//if (t1 == 780)
		//		//{
		//		//	cout << " CoordTest " << para->getParH(0)->coordX_SP[t1] << " " << para->getParH(0)->coordY_SP[t1] << " " << para->getParH(0)->coordZ_SP[t1] << " " << t1 << endl;
		//		//}
		//		//cout << " CoordTest " << para->getParH(0)->coordX_SP[t1] << " " << para->getParH(0)->coordY_SP[t1] << " " << para->getParH(0)->coordZ_SP[t1] << endl;
		//	}
		//}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//test2
		//cout << " erster Punkt des Vektors neighborX_SP: " << para->getParH(i)->neighborX_SP[1] << endl;
		//cout << " letzter Punkt des Vektors neighborX_SP: " << para->getParH(i)->neighborX_SP[para->getParH(i)->size_Mat_SP] << endl;
		//cout << " erster Punkt des Vektors neighborY_SP: " << para->getParH(i)->neighborY_SP[1] << endl;
		//cout << " letzter Punkt des Vektors neighborY_SP: " << para->getParH(i)->neighborY_SP[para->getParH(i)->size_Mat_SP] << endl;
		//cout << " erster Punkt des Vektors neighborZ_SP: " << para->getParH(i)->neighborZ_SP[1] << endl;
		//cout << " letzter Punkt des Vektors neighborZ_SP: " << para->getParH(i)->neighborZ_SP[para->getParH(i)->size_Mat_SP] << endl;
		//cout << " erster Punkt des Vektors geoSP: " << para->getParH(i)->geoSP[1] << endl;
		//cout << " letzter Punkt des Vektors geoSP: " << para->getParH(i)->geoSP[para->getParH(i)->size_Mat_SP] << endl;
		//////////////////////////////////////////////////////////////////////////
		for (int j = 0; j <= temp; j++)
		{
			para->getParH(i)->vx_SP[j] = para->getVelocity();//0.0f;//0.035f;
			para->getParH(i)->vy_SP[j] = 0.0f;//para->getVelocity();//0.0f;
			para->getParH(i)->vz_SP[j] = 0.0f;
			para->getParH(i)->rho_SP[j] = 0.0f;
			para->getParH(i)->press_SP[j] = 0.0f;
			if (para->getCalcMedian()){
				para->getParH(i)->vx_SP_Med[j] = 0.0f;
				para->getParH(i)->vy_SP_Med[j] = 0.0f;
				para->getParH(i)->vz_SP_Med[j] = 0.0f;
				para->getParH(i)->rho_SP_Med[j] = 0.0f;
				para->getParH(i)->press_SP_Med[j] = 0.0f;
			}
			if (para->getUseWale()) {
				para->getParH(i)->turbViscosity[j] = 0.0f;
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////Acoustic Test
		//doubflo factor01 = 20.0;
		//doubflo factor02 = 40.0;
		//doubflo factor03 = 60.0;
		//doubflo factor04 = 80.0;
		//doubflo factor = factor04; //01 = 20.0f; 02 = 40.0f; 03 = 60.0f; 04 = 80.0f

		//doubflo myEpsilon = 0.001; // * 0.1
		//doubflo b = 0.1 * factor; 
		//doubflo myAlpha = log(2.0) / (b * b);
		//doubflo myRadius = 0.0;
		//doubflo coordX = factor + 1.0;
		//doubflo coordY = factor + 1.0;
		//doubflo coordZ = factor + 1.0;
		//doubflo distX = 0.0;
		//doubflo distY = 0.0;
		//doubflo distZ = 0.0;
		//doubflo myPress = 0.0;

		//for (int j = 0; j <= temp; j++)
		//{
		//	distX = para->getParH(i)->coordX_SP[j] - coordX;
		//	distY = para->getParH(i)->coordY_SP[j] - coordY;
		//	distZ = para->getParH(i)->coordZ_SP[j] - coordZ;
		//	//myRadius = sqrt((distX * distX) + (distY * distY)); //XY
		//	//myRadius = sqrt((distY * distY) + (distZ * distZ)); //YZ
		//	myRadius = sqrt((distX * distX) + (distZ * distZ)); //XZ

		//	myPress = myEpsilon * exp(-myAlpha * myRadius * myRadius);

		//	para->getParH(i)->rho_SP[j] = myPress;
		//	para->getParH(i)->press_SP[j] = myPress;
		//}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Taylor Green Vortex uniform
		doubflo PI = 3.141592653589793238462643383279f;
		doubflo gridX = para->getParH(i)->gridNX - 1;
		doubflo gridY = para->getParH(i)->gridNY - 1;
		doubflo gridZ = para->getParH(i)->gridNZ - 1;
		//like MG
		//doubflo uAdvect = doubflo (1. / 250.); //32 nodes -> 250; 40 nodes -> 200; 64 nodes -> 500; 128 nodes -> 1000; 256 nodes -> 2000; 512 nodes -> 4000
		doubflo uAdvect = doubflo(0.0008); //32 nodes -> 0.032; 64 nodes -> 0.016; 128 nodes -> 0.008; 256 nodes -> 0.004; 512 nodes -> 0.002

		for (int j = 0; j <= temp; j++)
		{
			para->getParH(i)->rho_SP[j] = 
				(doubflo)((para->getVelocity()*para->getVelocity())*3.0 / 4.0*(cos(para->getParH(i)->coordX_SP[j]*4.0*PI / 
				(doubflo)gridX) + cos(para->getParH(i)->coordZ_SP[j] *4.0*PI /
				(doubflo)gridZ)))*(doubflo)(gridZ) / (doubflo)(gridX);
			
			para->getParH(i)->vy_SP[j] = (doubflo)0.0;

			//incl. additional velocity
			//para->getParH(i)->vx_SP[j] = 
			//	(doubflo)((32. * 32. ) / (1000. * 32.) * para->getVelocity() / 0.1 + para->getVelocity()*sin((para->getParH(i)->coordX_SP[j]*2.0*PI /
			//	(doubflo)gridX))*cos(para->getParH(i)->coordZ_SP[j]*2.0*PI / (doubflo)gridZ));
			//like MG
			para->getParH(i)->vx_SP[j] =
				(doubflo)(para->getVelocity()*sin((para->getParH(i)->coordX_SP[j] * 2.0*PI / (doubflo)gridX))*cos(para->getParH(i)->coordZ_SP[j] * 2.0*PI / (doubflo)gridZ)) + uAdvect * (1.0 + para->getParH(i)->rho_SP[j]);

			//no additional velocity
			//para->getParH(i)->vz_SP[j] = 
			//	(doubflo)(-para->getVelocity()*cos((para->getParH(i)->coordX_SP[j] *2.0*PI / (doubflo)gridX))*sin(para->getParH(i)->coordZ_SP[j]*2.0*PI /
			//	(doubflo)gridZ))*(doubflo)(gridZ) / (doubflo)(gridX);
			//like MG
			para->getParH(i)->vz_SP[j] =
				(doubflo)(-para->getVelocity()*cos((para->getParH(i)->coordX_SP[j] * 2.0*PI / (doubflo)gridX))*sin(para->getParH(i)->coordZ_SP[j] * 2.0*PI / (doubflo)gridZ));// *(doubflo)(gridZ) / (doubflo)(gridX);
			//cout << "\n Ich lebe: " << j << endl;
		}
		//cout << "\n GridNx: " << gridX << " GridNy: " << gridY << " GridNz: " << gridZ << endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////Desox fix
		//int levelToChange = 0;//uniform ... otherwise 7;
		//if (i==levelToChange)
		//{
		//	doubflo startX    = 769.0;//1025.0f;//uniform ... otherwise 114.999f; //120.95664f;
		//	doubflo endX      = 1537.0;//1793.0f;//uniform ... otherwise 120.999f;
		//	doubflo sizePress = endX - startX;
		//	doubflo tempPress = (1.0f / para->getFactorPressBC());
		//	cout << " tempPress: " << tempPress << endl;
		//	doubflo factor    = 0.0f;
		//	int tempIn = 0;
		//	for (int j = 0; j <= temp; j++)
		//	{
		//		if ( (para->getParH(i)->coordX_SP[j] > startX) && (para->getParH(i)->coordX_SP[j] < endX))
		//		{
		//			factor = (((para->getParH(i)->coordX_SP[j] - startX) / sizePress)); 
		//			//if (tempIn < 20)
		//			//{
		//				//cout << " CoordX: " << para->getParH(i)->coordX_SP[j] << endl;
		//				//cout << " factor: " << factor << endl;
		//			//}
		//			//tempIn++;
		//			para->getParH(i)->rho_SP[j] = factor * tempPress;
		//			para->getParH(i)->press_SP[j] = factor * tempPress;
		//		}
		//	}
		//}
		////////////////////////////////////////////////////////////////////////////////
		//Quatsch Test
		//for (int t1 = 0; t1 < 10; t1++)
		//{
		//	cout <<" CoordTest " << para->getParH(i)->coordX_SP[t1] << " " << para->getParH(i)->coordY_SP[t1] << " " << para->getParH(i)->coordZ_SP[t1] << endl;
		//}
		////////////////////////////////////////////////////////////////////////////
		rearrangeGeometry(para, i);
		////////////////////////////////////////////////////////////////////////////
		para->cudaCopySP(i);
		para->cudaCopyCoord(i);
		if (para->getCalcMedian()) 	    para->cudaCopyMedianSP(i);
		if (para->getCalcParticle()) 	para->cudaCopyNeighborWSB(i);
		if (para->getUseWale()) {
			para->cudaCopyNeighborWSB(i);
			para->cudaCopyTurbulentViscosityHD(i);
		}
		//cout <<"Test 6 " <<endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	cout << "Gesamtanzahl Knoten = " << AnzahlKnotenGes << endl;



	//////////////////////////////////////////////////////////////////////////////
	////Mist
	////////////////////////////////////////////////////////////////////////////
	//double minX, maxX;  
	//double minY, maxY;
	//double minZ, maxZ;
	//double diffX = 0.;
	//double diffY = 0.;
	//double diffZ = 0.;
	////int levelX, levelY, levelZ;

	//for (int i = 0; i < para->getMinCoordX().size(); i++)
	//{
	//	if ((para->getMaxCoordX()[i] - para->getMinCoordX()[i]) > diffX)
	//	{
	//		diffX = para->getMaxCoordX()[i] - para->getMinCoordX()[i];
	//		minX = para->getMinCoordX()[i];
	//		maxX = para->getMaxCoordX()[i];
	//		//levelX = i;
	//	} 
	//	if ((para->getMaxCoordY()[i] - para->getMinCoordY()[i]) > diffY)
	//	{
	//		diffY = para->getMaxCoordY()[i] - para->getMinCoordY()[i];
	//		minY = para->getMinCoordY()[i];
	//		maxY = para->getMaxCoordY()[i];
	//		//levelY = i;
	//	}
	//	if ((para->getMaxCoordZ()[i] - para->getMinCoordZ()[i]) > diffZ)
	//	{
	//		diffZ = para->getMaxCoordZ()[i] - para->getMinCoordZ()[i];
	//		minZ = para->getMinCoordZ()[i];
	//		maxZ = para->getMaxCoordZ()[i];
	//		//levelZ = i;
	//	}
	//}
	//
	//double midX = (maxX + minX) / 2.;
	//double midY = (maxY + minY) / 2.;
	//double midZ = (maxZ + minZ) / 2. + 10.0;
	////cout << " minX = " << minX << " minY = " << minY << " minZ = " << minZ << endl;
	////cout << " maxX = " << maxX << " maxY = " << maxY << " maxZ = " << maxZ << endl;
	////cout << " midx = " << midX << " midy = " << midY << " midz = " << midZ << endl;
	//double dx = 1.0;//(1.0/pow(2.0,i));

	////set velocity and press
	//for (int i=0; i<=level;i++) {
	//	int temp = coordX->getSize(i) + 1;
	//	//double dx = (1.0/pow(2.0,i));
	//	double d = maxY - minY + dx;
	//	double r = d / 2.0;
	//	double r2 = r * r;
	//	double l = (maxX - minX + dx);

	//	for (int j = 0; j <= temp; j++)
	//	{
	//		////////////////////////////////////////////////////////////////////////////
	//		para->getParH(i)->vz_SP[j]    = sin(3.14*2.0/diffX*para->getParH(i)->coordX_SP[j])*sin(3.14*2.0/diffY*para->getParH(i)->coordY_SP[j])*sin(3.14*2.0/diffZ*para->getParH(i)->coordZ_SP[j])*(0.000005f) / (3.0 * 4.0 * para->getViscosity() * l) * 
	//									    (r2 - ( ((midY - para->getParH(i)->coordY_SP[j]) * (midY - para->getParH(i)->coordY_SP[j])) + 
	//									    ((midZ - para->getParH(i)->coordZ_SP[j]) * (midZ - para->getParH(i)->coordZ_SP[j])) ));
	//		////////////////////////////////////////////////////////////////////////////
	//		//para->getParH(i)->rho_SP[j]   = (1.0 - ((para->getParH(i)->coordX_SP[j] - minX) / l * 2.0)) ;
	//		//////////////////////////////////////////////////////////////////////////////
	//		//para->getParH(i)->press_SP[j] = (1.0 - ((para->getParH(i)->coordX_SP[j] - minX) / l * 2.0)) ;
	//		////////////////////////////////////////////////////////////////////////////
	//	}
	//	////////////////////////////////////////////////////////////////////////////
	//	para->cudaCopySP(i);
	//	////////////////////////////////////////////////////////////////////////////
	//}
	////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	////Poiseuille Flow
	////////////////////////////////////////////////////////////////////////////
	//double minX, maxX;  
	//double minY, maxY;
	//double minZ, maxZ;
	//double diffX = 0.;
	//double diffY = 0.;
	//double diffZ = 0.;
	////int levelX, levelY, levelZ;

	//for (int i = 0; i < para->getMinCoordX().size(); i++)
	//{
	//	if ((para->getMaxCoordX()[i] - para->getMinCoordX()[i]) > diffX)
	//	{
	//		diffX = para->getMaxCoordX()[i] - para->getMinCoordX()[i];
	//		minX = para->getMinCoordX()[i];
	//		maxX = para->getMaxCoordX()[i];
	//		//levelX = i;
	//	} 
	//	if ((para->getMaxCoordY()[i] - para->getMinCoordY()[i]) > diffY)
	//	{
	//		diffY = para->getMaxCoordY()[i] - para->getMinCoordY()[i];
	//		minY = para->getMinCoordY()[i];
	//		maxY = para->getMaxCoordY()[i];
	//		//levelY = i;
	//	}
	//	if ((para->getMaxCoordZ()[i] - para->getMinCoordZ()[i]) > diffZ)
	//	{
	//		diffZ = para->getMaxCoordZ()[i] - para->getMinCoordZ()[i];
	//		minZ = para->getMinCoordZ()[i];
	//		maxZ = para->getMaxCoordZ()[i];
	//		//levelZ = i;
	//	}
	//}
	//
	//double midX = (maxX + minX) / 2.;
	//double midY = (maxY + minY) / 2.;
	//double midZ = (maxZ + minZ) / 2.;
	////cout << " minX = " << minX << " minY = " << minY << " minZ = " << minZ << endl;
	////cout << " maxX = " << maxX << " maxY = " << maxY << " maxZ = " << maxZ << endl;
	////cout << " midx = " << midX << " midy = " << midY << " midz = " << midZ << endl;
	//double dx = 1.0;//(1.0/pow(2.0,i));

	////set velocity and press
	//for (int i=0; i<=level;i++) {
	//	int temp = coordX->getSize(i) + 1;
	//	//double dx = (1.0/pow(2.0,i));
	//	double d = maxY - minY + dx;
	//	double r = d / 2.0;
	//	double r2 = r * r;
	//	double l = (maxX - minX + dx);

	//	for (int j = 0; j <= temp; j++)
	//	{
	//		////////////////////////////////////////////////////////////////////////////
	//		para->getParH(i)->vx_SP[j]    = (2.0 / para->getFactorPressBC()) / (3.0 * 4.0 * para->getViscosity() * l) * 
	//									    (r2 - ( ((midY - para->getParH(i)->coordY_SP[j]) * (midY - para->getParH(i)->coordY_SP[j])) + 
	//									    ((midZ - para->getParH(i)->coordZ_SP[j]) * (midZ - para->getParH(i)->coordZ_SP[j])) ));
	//		////////////////////////////////////////////////////////////////////////////
	//		para->getParH(i)->rho_SP[j]   = (1.0 - ((para->getParH(i)->coordX_SP[j] - minX) / l * 2.0)) / para->getFactorPressBC();
	//		////////////////////////////////////////////////////////////////////////////
	//		para->getParH(i)->press_SP[j] = (1.0 - ((para->getParH(i)->coordX_SP[j] - minX) / l * 2.0)) / para->getFactorPressBC();
	//		////////////////////////////////////////////////////////////////////////////
	//	}
	//	////////////////////////////////////////////////////////////////////////////
	//	para->cudaCopySP(i);
	//	////////////////////////////////////////////////////////////////////////////
	//}
	////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////
	//Coordinates
	//cout <<" X-Drehkoordinate: " << para->getParH(2)->coordX_SP[2822040] << endl;
	//cout <<" Y-Drehkoordinate: " << para->getParH(2)->coordY_SP[2822040] << endl;
	////////////////////////////////////////////////////////////////////////////
	////MS added
	//ifstream file_localToGlobal;
	//file_localToGlobal.open("D:/temp/Soeren/mesh_to_read/sphere_refinement2/coordLocalToGlobal.dat", ios::in);
	//if (!file_localToGlobal) {
	//	cerr << "Fehler beim Oeffnen von file_localToGlobal" <<endl;
	//	exit(1);
	//}
	////Startpositionen einlesen
	//string buffer;
	//for (int i=0; i<=level;i++) {
	//	file_localToGlobal > > buffer;
	//	file_localToGlobal > > para->getParH(i)->distX;
	//	file_localToGlobal > > para->getParH(i)->distY;
	//	file_localToGlobal > > para->getParH(i)->distZ;
	//}
	////////////////////////////////////////////////////////////////////////////

	//Test Soeren------------------------------------------------------------------------------------------------------------------------
	//cout << "Test Koordinanten: " << endl;
	//cout << "TEST coordX level: "<< i <<" , erster: " << ara_coordX[i][1] <<", letzter: "<< ara_coordX[i][temp-1] << endl;
	//cout << "TEST coordY level: "<< i <<" , erster: " << ara_coordY[i][1] <<", letzter: "<< ara_coordY[i][temp-1] << endl;
	//cout << "TEST coordZ level: "<< i <<" , erster: " << ara_coordZ[i][1] <<", letzter: "<< ara_coordZ[i][temp-1] << endl;
	//cout << "Test Nachbarn: " << endl;
	//cout << "TEST neighX level: "<< i <<" , erster: " << ara_nX[i][1] <<", letzter: "<< ara_nX[i][temp-1] << endl;
	//cout << "TEST neighY level: "<< i <<" , erster: " << ara_nY[i][1] <<", letzter: "<< ara_nY[i][temp-1] << endl;
	//cout << "TEST neighZ level: "<< i <<" , erster: " << ara_nZ[i][1] <<", letzter: "<< ara_nZ[i][temp-1] << endl;
	//cout << "Test geo: " << endl;
	//cout << "TEST geo level: "<< i <<" , erster: " << ara_geoV[i][1] <<", letzter: "<< ara_geoV[i][temp-1] << endl;
	//ENDE Test Soeren------------------------------------------------------------------------------------------------------------------------

	//}	

	//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------


	//Test Soeren------------------------------------------------------------------------------------------------------------------------
	//if(TEST==true) {
	//	cout << "Test Koordinanten: " << endl;
	//	cout << "TEST coordX level: "<< i <<" , erster: " << ara_coordX[i][1] <<", letzter: "<< ara_coordX[i][temp-1] << endl;
	//	cout << "TEST coordY level: "<< i <<" , erster: " << ara_coordY[i][1] <<", letzter: "<< ara_coordY[i][temp-1] << endl;
	//	cout << "TEST coordZ level: "<< i <<" , erster: " << ara_coordZ[i][1] <<", letzter: "<< ara_coordZ[i][temp-1] << endl;

	//	cout << "Test geo: " << endl;
	//	cout << "TEST geo level: "<< i <<" , erster: " << ara_geoV[i][1] <<", letzter: "<< ara_geoV[i][temp-1] << endl;
	////ENDE Test Soeren------------------------------------------------------------------------------------------------------------------------
	//}
	


	delete coordX;
	delete coordY;
	delete coordZ;

	delete geoV;

	cout << "-----Ende Coord, Neighbor, Geo------" <<endl;

}

void Interface::allocArrays_OffsetScale(Parameter* para) {
	
	cout << "-----Config Arrays OffsetScale------" <<endl;
	OffsetScale *obj_offCF = new OffsetScale(para->getscaleOffsetCF(), true);
	OffsetScale *obj_offFC = new OffsetScale(para->getscaleOffsetFC(), true);
	OffsetScale *obj_scaleCFC = new OffsetScale(para->getscaleCFC(), false);
	OffsetScale *obj_scaleCFF = new OffsetScale(para->getscaleCFF(), false);
	OffsetScale *obj_scaleFCC = new OffsetScale(para->getscaleFCC(), false);
	OffsetScale *obj_scaleFCF = new OffsetScale(para->getscaleFCF(), false);

	int level = obj_offCF->getLevel();

	int AnzahlKnotenGesCF = 0;
	int AnzahlKnotenGesFC = 0;

	for (int i=0; i<level; i++) {
		unsigned int tempCF = obj_offCF->getSize(i);
		cout << "Groesse der Daten CF vom Level "<< i <<" : " << tempCF <<endl;
		unsigned int tempFC = obj_offFC->getSize(i);
		cout << "Groesse der Daten FC vom Level "<< i <<" : " << tempFC <<endl;
		
		AnzahlKnotenGesCF += tempCF;
		AnzahlKnotenGesFC += tempFC;
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//size + memsize CF
		para->getParH(i)->K_CF = tempCF;
		para->getParD(i)->K_CF = para->getParH(i)->K_CF;
		para->getParH(i)->intCF.kCF = para->getParH(i)->K_CF;
		para->getParD(i)->intCF.kCF = para->getParH(i)->K_CF;
		para->getParH(i)->mem_size_kCF = sizeof(unsigned int)* para->getParH(i)->K_CF;
		para->getParD(i)->mem_size_kCF = sizeof(unsigned int)* para->getParD(i)->K_CF;
		para->getParH(i)->mem_size_kCF_off = sizeof(doubflo)* para->getParH(i)->K_CF;
		para->getParD(i)->mem_size_kCF_off = sizeof(doubflo)* para->getParD(i)->K_CF;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//size + memsize FC
		para->getParH(i)->K_FC = tempFC;
		para->getParD(i)->K_FC = para->getParH(i)->K_FC;
		para->getParH(i)->intFC.kFC = para->getParH(i)->K_FC;
		para->getParD(i)->intFC.kFC = para->getParH(i)->K_FC;
		para->getParH(i)->mem_size_kFC = sizeof(unsigned int)* para->getParH(i)->K_FC;
		para->getParD(i)->mem_size_kFC = sizeof(unsigned int)* para->getParD(i)->K_FC;
		para->getParH(i)->mem_size_kFC_off = sizeof(doubflo)* para->getParH(i)->K_FC;
		para->getParD(i)->mem_size_kFC_off = sizeof(doubflo)* para->getParD(i)->K_FC;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//alloc
		para->cudaAllocInterfaceCF(i);
		para->cudaAllocInterfaceFC(i);
		para->cudaAllocInterfaceOffCF(i);
		para->cudaAllocInterfaceOffFC(i);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//init
		obj_offCF->initArrayOffset(para->getParH(i)->offCF.xOffCF, para->getParH(i)->offCF.yOffCF, para->getParH(i)->offCF.zOffCF, i);
		obj_offFC->initArrayOffset(para->getParH(i)->offFC.xOffFC, para->getParH(i)->offFC.yOffFC, para->getParH(i)->offFC.zOffFC, i);
		obj_scaleCFC->initArray(para->getParH(i)->intCF.ICellCFC, i);
		obj_scaleCFF->initArray(para->getParH(i)->intCF.ICellCFF, i);
		obj_scaleFCC->initArray(para->getParH(i)->intFC.ICellFCC, i);
		obj_scaleFCF->initArray(para->getParH(i)->intFC.ICellFCF, i);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//copy
		para->cudaCopyInterfaceCF(i);
		para->cudaCopyInterfaceFC(i);
		para->cudaCopyInterfaceOffCF(i);
		para->cudaCopyInterfaceOffFC(i);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	cout << "Gesamtanzahl Knoten CF = " << AnzahlKnotenGesCF << endl;
	cout << "Gesamtanzahl Knoten FC = " << AnzahlKnotenGesFC << endl;
		
	delete obj_offCF;
	delete obj_offFC;
	delete obj_scaleCFC;
	delete obj_scaleCFF;
	delete obj_scaleFCC;
	delete obj_scaleFCF;
	cout << "-----Ende OffsetScale------" <<endl;
}

void Interface::allocArrays_BoundaryValues(Parameter* para) {

	cout << "-----Config Arrays BoundaryValues------" <<endl;
	int t = 6;


	BoundaryValues **BC_Values = new BoundaryValues*[t];
	this->sortSystem(BC_Values , para, t);

	way = new string[t];

	for (int i = 0; i < t; i++){
		this->way[i] = BC_Values[i]->getWay();
		cout << this->system[i] <<" Boundary: " << way[i] << endl;
	}

	BoundaryValues *obj_geomV=new BoundaryValues(para->getgeomBoundaryBcValues(), para, "geo");

	int level = BC_Values[0]->getLevel();


	vector<vector<vector<doubflo> > > pressureV;
	pressureV.resize(level+1);
	vector<vector<vector<doubflo> > > velocityV;
	velocityV.resize(level+1);
	//vector<vector<vector<unsigned int> > > periodV;
	//periodV.resize(level + 1);
	//vector<vector<unsigned int> > periodIndex;
	//periodIndex.resize(level + 1);
	vector<vector<vector<doubflo> > > outflowV;
	outflowV.resize(level+1);

	for(int i=0; i<=level; i++) {
		pressureV[i].resize(2); 
		velocityV[i].resize(3);
		//periodV[i].resize(1);
		outflowV[i].resize(2); 
	}

	//////////////////////////////////////////////////////////////////////////
	//added by Martin S.
	//////////////////////////////////////////////////////////////////////////
	////1D domain decomposition
	//vector< BoundaryValues* > procNeighborsSend;
	//vector< BoundaryValues* > procNeighborsRecv;
	//vector< int > neighborRank;

	//if (para->getNumprocs()>1)
	//{
	//	for (int i = 0; i < para->getNumprocs(); i++)
	//	{
	//		BoundaryValues *pnSend = new BoundaryValues(i, para, "send");
	//		BoundaryValues *pnRecv = new BoundaryValues(i, para, "recv");
	//		if (para->getIsNeighbor())
	//		{
	//			procNeighborsSend.push_back(pnSend);
	//			procNeighborsRecv.push_back(pnRecv);
	//			neighborRank.push_back(i);
	//		}
	//	}
	//}
	////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	vector< BoundaryValues* > procNeighborsSendX, procNeighborsSendY, procNeighborsSendZ;
	vector< BoundaryValues* > procNeighborsRecvX, procNeighborsRecvY, procNeighborsRecvZ;
	vector< int >             neighborRankX,      neighborRankY,      neighborRankZ;

	procNeighborsSendX.resize(0);
	procNeighborsSendY.resize(0);
	procNeighborsSendZ.resize(0);
	procNeighborsRecvX.resize(0);
	procNeighborsRecvY.resize(0);
	procNeighborsRecvZ.resize(0);
	neighborRankX.resize(0);
	neighborRankY.resize(0);
	neighborRankZ.resize(0);

	if (para->getNumprocs()>1)
	{
		for (int i = 0; i < para->getNumprocs(); i++)
		{
			BoundaryValues *pnXsend = new BoundaryValues(i, para, "send", "X");
			BoundaryValues *pnYsend = new BoundaryValues(i, para, "send", "Y");
			BoundaryValues *pnZsend = new BoundaryValues(i, para, "send", "Z");
			BoundaryValues *pnXrecv = new BoundaryValues(i, para, "recv", "X");
			BoundaryValues *pnYrecv = new BoundaryValues(i, para, "recv", "Y");
			BoundaryValues *pnZrecv = new BoundaryValues(i, para, "recv", "Z");
			if (para->getIsNeighborX())
			{
				procNeighborsSendX.push_back(pnXsend);
				procNeighborsRecvX.push_back(pnXrecv);
				neighborRankX.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankX: " << i << endl; 
			}
			if (para->getIsNeighborY())
			{
				procNeighborsSendY.push_back(pnYsend);
				procNeighborsRecvY.push_back(pnYrecv);
				neighborRankY.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankY: " << i << endl; 
			}
			if (para->getIsNeighborZ())
			{
				procNeighborsSendZ.push_back(pnZsend);
				procNeighborsRecvZ.push_back(pnZrecv);
				neighborRankZ.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankZ: " << i << endl; 
			}
		}
		cout << "MyID: " << para->getMyID() << ", size of neighborRankX: " << neighborRankX.size() << ", size of neighborRankY: " << neighborRankY.size() << ", size of neighborRankZ: " << neighborRankZ.size() << endl; 
	}
	////////////////////////////////////////////////////////////////////////
	//3D domain decomposition convection diffusion
	vector< BoundaryValues* > procNeighborsSendADX, procNeighborsSendADY, procNeighborsSendADZ;
	vector< BoundaryValues* > procNeighborsRecvADX, procNeighborsRecvADY, procNeighborsRecvADZ;
	vector< int >             neighborRankADX,      neighborRankADY,      neighborRankADZ;

	procNeighborsSendADX.resize(0);
	procNeighborsSendADY.resize(0);
	procNeighborsSendADZ.resize(0);
	procNeighborsRecvADX.resize(0);
	procNeighborsRecvADY.resize(0);
	procNeighborsRecvADZ.resize(0);
	neighborRankADX.resize(0);
	neighborRankADY.resize(0);
	neighborRankADZ.resize(0);

	if (para->getDiffOn()==true && para->getNumprocs()>1)
	{
		for (int i = 0; i < para->getNumprocs(); i++)
		{
			BoundaryValues *pnADXsend = new BoundaryValues(i, para, "send", "X");
			BoundaryValues *pnADYsend = new BoundaryValues(i, para, "send", "Y");
			BoundaryValues *pnADZsend = new BoundaryValues(i, para, "send", "Z");
			BoundaryValues *pnADXrecv = new BoundaryValues(i, para, "recv", "X");
			BoundaryValues *pnADYrecv = new BoundaryValues(i, para, "recv", "Y");
			BoundaryValues *pnADZrecv = new BoundaryValues(i, para, "recv", "Z");
			if (para->getIsNeighborX())
			{
				procNeighborsSendADX.push_back(pnADXsend);
				procNeighborsRecvADX.push_back(pnADXrecv);
				neighborRankADX.push_back(i);
			}
			if (para->getIsNeighborY())
			{
				procNeighborsSendADY.push_back(pnADYsend);
				procNeighborsRecvADY.push_back(pnADYrecv);
				neighborRankADY.push_back(i);
			}
			if (para->getIsNeighborZ())
			{
				procNeighborsSendADZ.push_back(pnADZsend);
				procNeighborsRecvADZ.push_back(pnADZrecv);
				neighborRankADZ.push_back(i);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//cout << "Laenge 1 periodV: " << periodV.size() << endl;
	//cout << "Laenge 2 periodV: " << periodV[0].size() << endl;
	//cout << "Laenge 3 periodV: " << periodV[0][0].size() << endl;


	////----------------------------------------Vektoren Initialisierung--------------
	string way_temp;
	for (int i = 0; i < t; i++){
		if (this->way[i] == "velocity") { velocityV = BC_Values[i]->setgetBoundarys(velocityV); }
		else if (this->way[i] == "pressure") { pressureV = BC_Values[i]->setgetBoundarys(pressureV); }
		else if (this->way[i] == "outflow") { outflowV = BC_Values[i]->setgetBoundarys(outflowV); }
		//else if (this->way[i] == "periodic_x" || this->way[i] == "periodic_y" || this->way[i] == "periodic_z") { periodV = BC_Values[i]->setgetBoundarys(periodV); way_temp = BC_Values[i]->getWay(); periodIndex = BC_Values[i]->setgetIndex(periodIndex); }

		//cout << "Laenge periodIndex: " << periodIndex[i].size() << endl;
		//cout << "Laenge 1 periodV: " << periodV.size() << endl;
		//cout << "Laenge 2 periodV: " << periodV[i].size() << endl;
		//cout << "Laenge 3 periodV: " << periodV[i][0].size() << endl;

		//if (way_temp != "")
		//{
		//	this->initPeriodicNeigh(periodV, periodIndex, way_temp);

		//		for (int i = 0; i<=level; i++)
		//		{
		//			////////////////////////////////////////////////////////////////////////
		//			neighX->initArray(para->getParH(i)->neighborX_SP,i);
		//			cout << " Test 3 " << endl;
		//			neighY->initArray(para->getParH(i)->neighborY_SP,i);
		//			cout << " Test 3 " << endl;
		//			neighZ->initArray(para->getParH(i)->neighborZ_SP,i);
		//			cout << " Test 3 " << endl;
		//			////////////////////////////////////////////////////////////////////////
		//			para->cudaCopySP(i);
		//			////////////////////////////////////////////////////////////////////////
		//		}

		//		cout << " Test 3 " << endl;

		//	//periodV wird für jede Seite resetet 
		//	for (int j = 0; j <= level; j++){
		//		periodV[j].clear();
		//	}
		//	periodV.clear();
		//	periodIndex.clear();

		//	cout << " Test 2 " << endl;
		//	periodV.resize(level + 1);
		//	periodIndex.resize(level + 1);
		//	for (int i = 0; i <= level; i++) {
		//		periodV[i].resize(1);
		//	}
		//	way_temp = "";

		//}
	}
	//cout << " Test 4 " << endl;

	////------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------
	//
	////------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------//------------------------



	////Soeren TEST---------------------------------------------------------------------------------------------------------------------------------------------------
	//if (TEST == true) {
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)velocityV[i][0].size();
	//		cout << "Groesse Velocity Values lvl " << i << " : " << temp << endl;
	//	}
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)pressureV[i][0].size();
	//		cout << "Groesse Pressure Values lvl " << i << " : " << temp << endl;
	//	}

	//	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	int ausgabe = 5;
	//	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	for (int j = 0; j <= level; j++) {
	//		if (pressureV[j][0].size() > 0){
	//			for (int i = 0; i < ausgabe; i++) {
	//				cout << "Pressure Values level " << j << ", eintrag " << i << " : " << pressureV[j][0][i] << " | " << pressureV[j][1][i] << endl;
	//			}
	//			for (unsigned int i = pressureV[j][0].size() - ausgabe; i < pressureV[j][0].size(); i++) {
	//				cout << "Pressure Values level " << j << ", eintrag " << i << " : " << pressureV[j][0][i] << " | " << pressureV[j][1][i] << endl;
	//			}
	//		}
	//	}
	//	for (int j = 0; j <= level; j++) {
	//		if (velocityV[j][0].size() >0){
	//			for (int i = 0; i < ausgabe; i++) {
	//				cout << "Velocity Values level " << j << ", eintrag " << i << " : " << velocityV[j][0][i] << " | " << velocityV[j][1][i] << " | " << velocityV[j][2][i] << endl;
	//			}
	//			for (unsigned int i = velocityV[j][0].size() - ausgabe; i < velocityV[j][0].size(); i++) {
	//				cout << "velocity Values level " << j << ", eintrag " << i << " : " << velocityV[j][0][i] << " | " << velocityV[j][1][i] << " | " << velocityV[j][2][i] << endl;
	//			}
	//		}
	//	}
	//	


	//	//TEST NACHBARN

	//	//for (int i = 300; i <= 500; i++) {
	//	//	cout << "nachbarX index " << i << " : " << neighX->getVec(0)[i] << endl;
	//	//}

	//	//for (int i = 300; i <= 800; i++) {
	//	//	cout << "nachbarY index " << i << " : " << neighY->getVec(0)[i] << endl;
	//	//}

	//}


	//ENDE TEST Soeren---------------------------------------------------------------------------------------------------------------------------------------------------
	



	//---------------------------------------------------------------------//AB HIER CODE FUER MARTIN////-------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp1 = (int)pressureV[i][0].size();
		if (temp1 > 1)
		{
			cout << "Groesse pressure level " << i << " : " << temp1 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->getParH(i)->QPress.kQ = temp1;
			para->getParD(i)->QPress.kQ = temp1;
			para->getParH(i)->kPressQread = temp1 * para->getD3Qxx();
			para->getParD(i)->kPressQread = temp1 * para->getD3Qxx();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocPress(i);
			///////////////////////////////
			////only for round of error test
			//para->cudaAllocTestRE(i, temp1);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int d = 0;
			int j = 0;
			int n = 0;

			for (vector<vector<vector<doubflo> > >::iterator it = pressureV.begin(); it != pressureV.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							if (j == 0) para->getParH(i)->QPress.RhoBC[n] = *it3;
							if (j == 1) para->getParH(i)->QPress.kN[n] = (int)*it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp1; m++)
			{
				para->getParH(i)->QPress.RhoBC[m] = (para->getParH(i)->QPress.RhoBC[m] / para->getFactorPressBC());
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyPress(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////only for round of error test
			//para->cudaCopyTestREtoDevice(i,temp1);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//// advection - diffusion stuff
			////cout << "vor advec diff" << endl;
			//if (para->getDiffOn()==true){
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor setzen von kTemp" << endl;
			//	para->getParH(i)->TempPress.kTemp = temp1;
			//	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	para->cudaAllocTempPressBC(i);
			//	//cout << "nach alloc" << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	for (int m = 0; m < temp1; m++)
			//	{
			//		para->getParH(i)->TempPress.temp[m] = para->getTemperatureInit();
			//		para->getParH(i)->TempPress.velo[m] = (doubflo)0.0;
			//		para->getParH(i)->TempPress.k[m]    = para->getParH(i)->QPress.k[m];
			//	}
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor copy" << endl;
			//	para->cudaCopyTempPressBCHD(i);
			//	//cout << "nach copy" << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife




	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp2 = (int)velocityV[i][0].size();
		if (temp2 > 1)
		{
			cout << "Groesse velocity level " << i << " : " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			int blocks = (temp2 / para->getParH(i)->numberofthreads) + 1;
			para->getParH(i)->Qinflow.kArray = blocks * para->getParH(i)->numberofthreads;
			para->getParD(i)->Qinflow.kArray = para->getParH(i)->Qinflow.kArray;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->getParH(i)->Qinflow.kQ = temp2;
			para->getParD(i)->Qinflow.kQ = temp2;
			para->getParH(i)->kInflowQ = temp2;
			para->getParD(i)->kInflowQ = temp2;
			para->getParH(i)->kInflowQread = temp2 * para->getD3Qxx();
			para->getParD(i)->kInflowQread = temp2 * para->getD3Qxx();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocVeloBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = velocityV.begin(); it != velocityV.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							if (j == 0) para->getParH(i)->Qinflow.Vx[n] = *it3;
							if (j == 1) para->getParH(i)->Qinflow.Vy[n] = *it3;
							if (j == 2) para->getParH(i)->Qinflow.Vz[n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						//cout << "n = " << n << endl;
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			//cout << "temp2 = " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp2; m++)
			{
				//para->getParH(i)->Qinflow.Vx[m] = para->getParH(i)->Qinflow.Vx[m] / para->getVelocityRatio();
				//para->getParH(i)->Qinflow.Vy[m] = para->getParH(i)->Qinflow.Vy[m] / para->getVelocityRatio();
				//para->getParH(i)->Qinflow.Vz[m] = para->getParH(i)->Qinflow.Vz[m] / para->getVelocityRatio();
				para->getParH(i)->Qinflow.Vx[m] = para->getVelocity();//0.035;
				para->getParH(i)->Qinflow.Vy[m] = 0.0;//para->getVelocity();//0.0;
				para->getParH(i)->Qinflow.Vz[m] = 0.0;
				//if (para->getParH(i)->Qinflow.Vz[m] > 0)
				//{
				//	cout << "velo Z = " << para->getParH(i)->Qinflow.Vz[m] << endl;
				//}
			}
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyVeloBC(i);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//// advection - diffusion stuff
			//if (para->getDiffOn()==true){
			//	//////////////////////////////////////////////////////////////////////////
			//	para->getParH(i)->TempVel.kTemp = temp2;
			//	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
			//	cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
			//	cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	para->cudaAllocTempVeloBC(i);
			//	//cout << "nach alloc " << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	for (int m = 0; m < temp2; m++)
			//	{
			//		para->getParH(i)->TempVel.temp[m]      = para->getTemperatureInit();
			//		para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
			//		para->getParH(i)->TempVel.velo[m]      = para->getVelocity();
			//		para->getParH(i)->TempVel.k[m]         = para->getParH(i)->Qinflow.k[m];
			//	}
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor copy " << endl;
			//	para->cudaCopyTempVeloBCHD(i);
			//	//cout << "nach copy " << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife


	//cout << "Test 1 " << endl;

	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp = (int)outflowV[i][0].size();
		if (temp > 1)
		{
			cout << "Groesse outflow level " << i << " : " << temp << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->getParH(i)->Qoutflow.kQ = temp;
			para->getParD(i)->Qoutflow.kQ = temp;
			para->getParH(i)->kOutflowQread = temp * para->getD3Qxx();
			para->getParD(i)->kOutflowQread = temp * para->getD3Qxx();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocOutflowBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int d = 0;
			int j = 0;
			int n = 0;

			for (vector<vector<vector<doubflo> > >::iterator it = outflowV.begin(); it != outflowV.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							if (j == 0) para->getParH(i)->Qoutflow.RhoBC[n] = *it3;
							if (j == 1) para->getParH(i)->Qoutflow.kN[n] = (int)*it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp; m++)
			{
				para->getParH(i)->Qoutflow.RhoBC[m] = (para->getParH(i)->Qoutflow.RhoBC[m] / para->getFactorPressBC()) * (doubflo)0.0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyOutflowBC(i);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife




	//--------------------------------------------------------------------------//
	if (para->getIsGeometryValues()){
		for (int i = 0; i <= level; i++) {
			int temp4 = obj_geomV->getSize(i);
			if (temp4 > 0)
			{
				cout << "Groesse der Daten obj_geomV, Level " << i << " : " << temp4 << endl;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->getParH(i)->QGeom.kQ = temp4;
				para->getParD(i)->QGeom.kQ = temp4;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaAllocGeomValuesBC(i);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_geomV->initArray(para->getParH(i)->QGeom.Vx, i, 0);
				obj_geomV->initArray(para->getParH(i)->QGeom.Vy, i, 1);
				obj_geomV->initArray(para->getParH(i)->QGeom.Vz, i, 2);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp4; m++)
				{
					para->getParH(i)->QGeom.Vx[m] = para->getParH(i)->QGeom.Vx[m] / para->getVelocityRatio();
					para->getParH(i)->QGeom.Vy[m] = para->getParH(i)->QGeom.Vy[m] / para->getVelocityRatio();
					para->getParH(i)->QGeom.Vz[m] = para->getParH(i)->QGeom.Vz[m] / para->getVelocityRatio();
					//para->getParH(i)->QGeom.Vx[m] = para->getParH(i)->QGeom.Vx[m] / 100.0f;
					//para->getParH(i)->QGeom.Vy[m] = para->getParH(i)->QGeom.Vy[m] / 100.0f;
					//para->getParH(i)->QGeom.Vz[m] = para->getParH(i)->QGeom.Vz[m] / 100.0f;
					//para->getParH(i)->QGeom.Vx[m] = 0.0f;
					//para->getParH(i)->QGeom.Vy[m] = 0.0f;
					//para->getParH(i)->QGeom.Vz[m] = 0.0f;
				}
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////Täst
				//for (int m = 0; m < temp4; m++)
				//{
				//	para->getParH(i)->QGeom.Vx[m] = para->getVelocity();//0.035f;
				//	para->getParH(i)->QGeom.Vy[m] = 0.0f;//para->getVelocity();//0.0f;
				//	para->getParH(i)->QGeom.Vz[m] = 0.0f;
				//}
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaCopyGeomValuesBC(i);
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//// advection - diffusion stuff
				//if (para->getDiffOn()==true){
				//	//////////////////////////////////////////////////////////////////////////
				//	para->getParH(i)->Temp.kTemp = temp4;
				//	cout << "Groesse kTemp = " << para->getParH(i)->Temp.kTemp << endl;
				//	//////////////////////////////////////////////////////////////////////////
				//	para->cudaAllocTempNoSlipBC(i);
				//	//////////////////////////////////////////////////////////////////////////
				//	for (int m = 0; m < temp4; m++)
				//	{
				//		para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
				//		para->getParH(i)->Temp.k[m]    = para->getParH(i)->QGeom.k[m];
				//	}
				//	//////////////////////////////////////////////////////////////////////////
				//	para->cudaCopyTempNoSlipBCHD(i);
				//	//////////////////////////////////////////////////////////////////////////
				//}
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
		}
	}//ende geo

	//cout << "Test 2 " << endl;

	////--------------------------------------------------------------------------//
	//if (para->getIsProp()){
	//	BoundaryValues *obj_propV=new BoundaryValues(para->getpropellerValues(), para, "prop");
	//	for (int i = 0; i <= level; i++) {
	//		int temp4 = obj_propV->getSize(i);
	//		if (temp4 > 0)
	//		{
	//			cout << "Groesse der Daten PropellerValues, Level " << i << " : " << temp4 << endl;
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			para->getParH(i)->QPropeller.kQ = temp4;
	//			para->getParD(i)->QPropeller.kQ = temp4;
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			para->cudaAllocVeloPropeller(i);
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			//Indexarray
	//			obj_propV->initIndex(para->getParH(i)->QPropeller.k, i);
	//			obj_propV->initArray(para->getParH(i)->QPropeller.Vx, i, 0);
	//			obj_propV->initArray(para->getParH(i)->QPropeller.Vy, i, 1);
	//			obj_propV->initArray(para->getParH(i)->QPropeller.Vz, i, 2);
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			for (int m = 0; m < temp4; m++)
	//			{
	//				para->getParH(i)->QPropeller.Vx[m] = para->getParH(i)->QPropeller.Vx[m] / para->getVelocityRatio();
	//				para->getParH(i)->QPropeller.Vy[m] = para->getParH(i)->QPropeller.Vy[m] / para->getVelocityRatio();
	//				para->getParH(i)->QPropeller.Vz[m] = para->getParH(i)->QPropeller.Vz[m] / para->getVelocityRatio();
	//			}
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			para->cudaCopyVeloPropeller(i);
	//			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		}
	//	}
	//}//ende prop

	//cout << "Test 3 " << endl;


	//--------------------------------------------------------------------------//
	//BoundaryValues *obj_cpTop=new BoundaryValues(para->getcpTop(), para, "cp");
	//BoundaryValues *obj_cpBottom=new BoundaryValues(para->getcpBottom(), para, "cp");
	//BoundaryValues *obj_cpBottom2=new BoundaryValues(para->getcpBottom2(), para, "cp");
	//if (para->getIsCp()){
	//	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	//Top
	//	for (int i = 0; i <= level; i++) {
	//		int temp = obj_cpTop->getSize(i);
	//		if (temp > 0)
	//		{
	//			cout << "Groesse der Daten CpTop, Level " << i << " : " << temp << endl;
	//			////////////////////////////////////////////////////////////////////////////
	//			para->getParH(i)->numberOfPointsCpTop = temp;
	//			para->getParD(i)->numberOfPointsCpTop = temp;
	//			////////////////////////////////////////////////////////////////////////////
	//			para->cudaAllocCpTop(i);
	//			////////////////////////////////////////////////////////////////////////////
	//			//Indexarray
	//			obj_cpTop->initIndex(para->getParH(i)->cpTopIndex, i);
	//			////////////////////////////////////////////////////////////////////////////
	//			for (int m = 0; m < temp; m++)
	//			{
	//				para->getParH(i)->cpPressTop[m] = 0.0;
	//			}
	//			////////////////////////////////////////////////////////////////////////////
	//			para->cudaCopyCpTopInit(i);
	//			////////////////////////////////////////////////////////////////////////////
	//		}
	//	}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////Bottom
		//for (int i = 0; i <= level; i++) {
		//	int temp = obj_cpBottom->getSize(i);
		//	if (temp > 0)
		//	{
		//		cout << "Groesse der Daten CpBottom, Level " << i << " : " << temp << endl;
		//		////////////////////////////////////////////////////////////////////////////
		//		para->getParH(i)->numberOfPointsCpBottom = temp;
		//		para->getParD(i)->numberOfPointsCpBottom = temp;
		//		////////////////////////////////////////////////////////////////////////////
		//		para->cudaAllocCpBottom(i);
		//		////////////////////////////////////////////////////////////////////////////
		//		//Indexarray
		//		obj_cpBottom->initIndex(para->getParH(i)->cpBottomIndex, i);
		//		////////////////////////////////////////////////////////////////////////////
		//		for (int m = 0; m < temp; m++)
		//		{
		//			para->getParH(i)->cpPressBottom[m] = 0.0;
		//		}
		//		////////////////////////////////////////////////////////////////////////////
		//		para->cudaCopyCpBottomInit(i);
		//		////////////////////////////////////////////////////////////////////////////
		//	}
		//}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////Bottom 2
		//for (int i = 0; i <= level; i++) {
		//	int temp = obj_cpBottom2->getSize(i);
		//	if (temp > 0)
		//	{
		//		cout << "Groesse der Daten CpBottom2, Level " << i << " : " << temp << endl;
		//		////////////////////////////////////////////////////////////////////////////
		//		para->getParH(i)->numberOfPointsCpBottom2 = temp;
		//		para->getParD(i)->numberOfPointsCpBottom2 = temp;
		//		////////////////////////////////////////////////////////////////////////////
		//		para->cudaAllocCpBottom2(i);
		//		////////////////////////////////////////////////////////////////////////////
		//		//Indexarray
		//		obj_cpBottom2->initIndex(para->getParH(i)->cpBottom2Index, i);
		//		////////////////////////////////////////////////////////////////////////////
		//		for (int m = 0; m < temp; m++)
		//		{
		//			para->getParH(i)->cpPressBottom2[m] = 0.0;
		//		}
		//		////////////////////////////////////////////////////////////////////////////
		//		para->cudaCopyCpBottom2Init(i);
		//		////////////////////////////////////////////////////////////////////////////
		//	}
		//}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	delete obj_cpTop;
	//	//delete obj_cpBottom;
	//	//delete obj_cpBottom2;
	//}//ende cp

	//cout << "Test 4 " << endl;


	//--------------------------------------------------------------------------//
	if (para->getConcFile()){
		BoundaryValues *obj_Conc=new BoundaryValues(para->getConcentration());
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//concentration
		for (int i = 0; i <= level; i++) {
			int temp = obj_Conc->getSize(i);
			if (temp > 0)
			{
				cout << "Groesse der Daten Concentration, Level " << i << " : " << temp << endl;
				////////////////////////////////////////////////////////////////////////////
				para->getParH(i)->numberOfPointsConc = temp;
				para->getParD(i)->numberOfPointsConc = temp;
				////////////////////////////////////////////////////////////////////////////
				para->cudaAllocConcFile(i);
				////////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_Conc->initIndex(para->getParH(i)->concIndex, i);
				////////////////////////////////////////////////////////////////////////////
				para->cudaCopyConcFile(i);
				////////////////////////////////////////////////////////////////////////////
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		delete obj_Conc;
	}//end concentration

	//cout << "Test 5 " << endl;



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//processor boundary (Martin Sch.) 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////1D domain decomposition
	//if ( (para->getNumprocs() > 1) && (procNeighborsSend.size() == procNeighborsRecv.size()) )
	//{
	//	for (int j = 0; j < procNeighborsSend.size(); j++)
	//	{
	//		for (int i = 0; i <= level; i++) {
	//			int tempSend = procNeighborsSend[j]->getSize(i);
	//			int tempRecv = procNeighborsRecv[j]->getSize(i);
	//			if (tempSend > 0)
	//			{
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//send
	//				cout << "Groesse der Daten für den Sendepuffer, Level " << i << " : " << tempSend << endl;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				para->setNumberOfProcessNeighbors((unsigned int)procNeighborsSend.size(), i, "send");
	//				para->getParH(i)->sendProcessNeighbor[j].rankNeighbor = neighborRank[j];
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				para->getParH(i)->sendProcessNeighbor[j].numberOfNodes = tempSend;
	//				para->getParD(i)->sendProcessNeighbor[j].numberOfNodes = tempSend;
	//				para->getParH(i)->sendProcessNeighbor[j].numberOfFs    = para->getD3Qxx() * tempSend;
	//				para->getParD(i)->sendProcessNeighbor[j].numberOfFs    = para->getD3Qxx() * tempSend;
	//				para->getParH(i)->sendProcessNeighbor[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
	//				para->getParD(i)->sendProcessNeighbor[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
	//				para->getParH(i)->sendProcessNeighbor[j].memsizeFs     = sizeof(doubflo)     *tempSend;
	//				para->getParD(i)->sendProcessNeighbor[j].memsizeFs     = sizeof(doubflo)     *tempSend;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//cout << "Groesse send numberOfNodes Host, Level " << i << " : " << para->getParH(i)->sendProcessNeighbor[j].numberOfNodes << endl;
	//				//cout << "Groesse send numberOfNodes Device, Level " << i << " : " << para->getParD(i)->sendProcessNeighbor[j].numberOfNodes << endl;
	//				//cout << "Groesse send numberOfFs Host, Level " << i << " : " << para->getParH(i)->sendProcessNeighbor[j].numberOfFs << endl;
	//				//cout << "Groesse send numberOfFs Device, Level " << i << " : " << para->getParD(i)->sendProcessNeighbor[j].numberOfFs << endl;
	//				//cout << "Groesse send memsizeIndex Host, Level " << i << " : " << para->getParH(i)->sendProcessNeighbor[j].memsizeIndex << endl;
	//				//cout << "Groesse send memsizeIndex Device, Level " << i << " : " << para->getParD(i)->sendProcessNeighbor[j].memsizeIndex << endl;
	//				//cout << "Groesse send memsizeFs Host, Level " << i << " : " << para->getParH(i)->sendProcessNeighbor[j].memsizeFs << endl;
	//				//cout << "Groesse send memsizeFs Device, Level " << i << " : " << para->getParD(i)->sendProcessNeighbor[j].memsizeFs << endl;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//recv
	//				cout << "Groesse der Daten für den Empfangspuffer, Level " << i << " : " << tempRecv << endl;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				para->setNumberOfProcessNeighbors((unsigned int)procNeighborsRecv.size(), i, "recv");
	//				para->getParH(i)->recvProcessNeighbor[j].rankNeighbor = neighborRank[j];
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				para->getParH(i)->recvProcessNeighbor[j].numberOfNodes = tempRecv;
	//				para->getParD(i)->recvProcessNeighbor[j].numberOfNodes = tempRecv;
	//				para->getParH(i)->recvProcessNeighbor[j].numberOfFs    = para->getD3Qxx() * tempRecv;
	//				para->getParD(i)->recvProcessNeighbor[j].numberOfFs    = para->getD3Qxx() * tempRecv;
	//				para->getParH(i)->recvProcessNeighbor[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
	//				para->getParD(i)->recvProcessNeighbor[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
	//				para->getParH(i)->recvProcessNeighbor[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
	//				para->getParD(i)->recvProcessNeighbor[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//cout << "Groesse recv numberOfNodes Host, Level " << i << " : " << para->getParH(i)->recvProcessNeighbor[j].numberOfNodes << endl;
	//				//cout << "Groesse recv numberOfNodes Device, Level " << i << " : " << para->getParD(i)->recvProcessNeighbor[j].numberOfNodes << endl;
	//				//cout << "Groesse recv numberOfFs Host, Level " << i << " : " << para->getParH(i)->recvProcessNeighbor[j].numberOfFs << endl;
	//				//cout << "Groesse recv numberOfFs Device, Level " << i << " : " << para->getParD(i)->recvProcessNeighbor[j].numberOfFs << endl;
	//				//cout << "Groesse recv memsizeIndex Host, Level " << i << " : " << para->getParH(i)->recvProcessNeighbor[j].memsizeIndex << endl;
	//				//cout << "Groesse recv memsizeIndex Device, Level " << i << " : " << para->getParD(i)->recvProcessNeighbor[j].memsizeIndex << endl;
	//				//cout << "Groesse recv memsizeFs Host, Level " << i << " : " << para->getParH(i)->recvProcessNeighbor[j].memsizeFs << endl;
	//				//cout << "Groesse recv memsizeFs Device, Level " << i << " : " << para->getParD(i)->recvProcessNeighbor[j].memsizeFs << endl;
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//malloc on host and device
	//				para->cudaAllocProcessNeighbor(i, j);
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				//init index arrays
	//				procNeighborsSend[j]->initIndex(para->getParH(i)->sendProcessNeighbor[j].index, i);
	//				procNeighborsRecv[j]->initIndex(para->getParH(i)->recvProcessNeighbor[j].index, i);
	//				////////////////////////////////////////////////////////////////////////////////////////
	//				para->cudaCopyProcessNeighborIndex(i, j);
	//				////////////////////////////////////////////////////////////////////////////////////////
	//			}
	//		}
	//	}
	//}//ende processor boundarys
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	// X
	if ( (para->getNumprocs() > 1) && (procNeighborsSendX.size() == procNeighborsRecvX.size()) )
	{
		for (int j = 0; j < procNeighborsSendX.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendX[j]->getSize(i);
				int tempRecv = procNeighborsRecvX[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den X Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsSendX.size(), i, "send");
					para->getParH(i)->sendProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborX[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den X Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsRecvX.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborX(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendX[j]->initIndex(para->getParH(i)->sendProcessNeighborX[j].index, i);
					procNeighborsRecvX[j]->initIndex(para->getParH(i)->recvProcessNeighborX[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborXIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende X processor boundarys
	//////////////////////////////////////////////////////////////////////////
	// Y
	if ( (para->getNumprocs() > 1) && (procNeighborsSendY.size() == procNeighborsRecvY.size()) )
	{
		for (int j = 0; j < procNeighborsSendY.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendY[j]->getSize(i);
				int tempRecv = procNeighborsRecvY[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den Y Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsSendY.size(), i, "send");
					para->getParH(i)->sendProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborY[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den Y Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsRecvY.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborY(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendY[j]->initIndex(para->getParH(i)->sendProcessNeighborY[j].index, i);
					procNeighborsRecvY[j]->initIndex(para->getParH(i)->recvProcessNeighborY[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborYIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende Y processor boundarys
	//////////////////////////////////////////////////////////////////////////
	// Z
	if ( (para->getNumprocs() > 1) && (procNeighborsSendZ.size() == procNeighborsRecvZ.size()) )
	{
		for (int j = 0; j < procNeighborsSendZ.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendZ[j]->getSize(i);
				int tempRecv = procNeighborsRecvZ[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den Z Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsSendZ.size(), i, "send");
					para->getParH(i)->sendProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfFs    = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeFs     = sizeof(doubflo)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den Z Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsRecvZ.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfFs    = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeFs    = sizeof(doubflo)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborZ(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendZ[j]->initIndex(para->getParH(i)->sendProcessNeighborZ[j].index, i);
					procNeighborsRecvZ[j]->initIndex(para->getParH(i)->recvProcessNeighborZ[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborZIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende Z processor boundarys
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//3D domain decomposition convection diffusion
	if (para->getDiffOn()==true){
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// X
		if ( (para->getNumprocs() > 1) && (procNeighborsSendADX.size() == procNeighborsRecvADX.size()) )
		{
			for (int j = 0; j < procNeighborsSendADX.size(); j++)
			{
				for (int i = 0; i <= level; i++) {
					int tempSend = procNeighborsSendADX[j]->getSize(i);
					int tempRecv = procNeighborsRecvADX[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADX[j].rankNeighbor = neighborRankADX[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADX[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADX[j].rankNeighbor = neighborRankADX[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADX[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADX(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADX[j]->initIndex(para->getParH(i)->sendProcessNeighborADX[j].index, i);
						procNeighborsRecvADX[j]->initIndex(para->getParH(i)->recvProcessNeighborADX[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADXIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende X processor boundarys
		//////////////////////////////////////////////////////////////////////////
		// Y
		if ( (para->getNumprocs() > 1) && (procNeighborsSendADY.size() == procNeighborsRecvADY.size()) )
		{
			for (int j = 0; j < procNeighborsSendADY.size(); j++)
			{
				for (int i = 0; i <= level; i++) {
					int tempSend = procNeighborsSendADY[j]->getSize(i);
					int tempRecv = procNeighborsRecvADY[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADY[j].rankNeighbor = neighborRankADY[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADY[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADY[j].rankNeighbor = neighborRankADY[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADY[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADY(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADY[j]->initIndex(para->getParH(i)->sendProcessNeighborADY[j].index, i);
						procNeighborsRecvADY[j]->initIndex(para->getParH(i)->recvProcessNeighborADY[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADYIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende Y processor boundarys
		//////////////////////////////////////////////////////////////////////////
		// Z
		if ( (para->getNumprocs() > 1) && (procNeighborsSendADZ.size() == procNeighborsRecvADZ.size()) )
		{
			for (int j = 0; j < procNeighborsSendADZ.size(); j++)
			{
				for (int i = 0; i <= level; i++) {
					int tempSend = procNeighborsSendADZ[j]->getSize(i);
					int tempRecv = procNeighborsRecvADZ[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADZ[j].rankNeighbor = neighborRankADZ[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADZ[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].numberOfFs    = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].memsizeIndex  = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].memsizeFs     = sizeof(doubflo)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADZ[j].rankNeighbor = neighborRankADZ[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADZ[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].numberOfFs    = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].memsizeIndex  = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].memsizeFs     = sizeof(doubflo)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADZ(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADZ[j]->initIndex(para->getParH(i)->sendProcessNeighborADZ[j].index, i);
						procNeighborsRecvADZ[j]->initIndex(para->getParH(i)->recvProcessNeighborADZ[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADZIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende Z processor boundarys
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}




	//---------------------------------------------------------------------//BIS HIER CODE FUER MARTIN////-------------------------------------------------------------------//
	

	for (int i = 0; i < t; i++) {
		delete BC_Values[i];
	}

	delete BC_Values;
	delete obj_geomV;

	cout << "-----Ende config Arrays BoundaryValues-----" <<endl;
}

void Interface::allocArrays_BoundaryQs(Parameter* para) {
	cout << "-----Config Arrays BoundaryQs------" <<endl;
	//cout << "1: MyID: " << para->getMyID() << endl;
	int t = 6;

	//cout << "2: MyID: " << para->getMyID() << endl;
	BoundaryQs **BC_Qs = new BoundaryQs*[t];
	this->sortSystem(BC_Qs, para, t);

	
	//cout << "3: MyID: " << para->getMyID() << endl;
	BoundaryQs *obj_geomQ = new BoundaryQs(para->getgeomBoundaryBcQs(), para, "geo", false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Normals Geo
	BoundaryQs *obj_geomNormalX = new BoundaryQs(para->getgeomBoundaryNormalX(), para, "geoNormal", false);
	BoundaryQs *obj_geomNormalY = new BoundaryQs(para->getgeomBoundaryNormalY(), para, "geoNormal", false);
	BoundaryQs *obj_geomNormalZ = new BoundaryQs(para->getgeomBoundaryNormalZ(), para, "geoNormal", false);
	//////////////////////////////////////////////////////////////////////////
	//Normals Inflow
	BoundaryQs *obj_inflowNormalX = new BoundaryQs(para->getInflowBoundaryNormalX(), para, "inflowNormal", false);
	BoundaryQs *obj_inflowNormalY = new BoundaryQs(para->getInflowBoundaryNormalY(), para, "inflowNormal", false);
	BoundaryQs *obj_inflowNormalZ = new BoundaryQs(para->getInflowBoundaryNormalZ(), para, "inflowNormal", false);
	//////////////////////////////////////////////////////////////////////////
	//Normals Outflow
	BoundaryQs *obj_outflowNormalX = new BoundaryQs(para->getOutflowBoundaryNormalX(), para, "outflowNormal", false);
	BoundaryQs *obj_outflowNormalY = new BoundaryQs(para->getOutflowBoundaryNormalY(), para, "outflowNormal", false);
	BoundaryQs *obj_outflowNormalZ = new BoundaryQs(para->getOutflowBoundaryNormalZ(), para, "outflowNormal", false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int level = BC_Qs[0]->getLevel();
	
	//cout << "4: MyID: " << para->getMyID() << endl;

	//-----------------------------------------Vektoren Deklarationen-----------
	vector<vector<vector<doubflo> > > noslipQs;
	vector<vector<unsigned int> > noslipIndex;
	vector<vector<vector<doubflo> > > slipQs;
	vector<vector<unsigned int> > slipIndex;
	vector<vector<vector<doubflo> > > pressureQs;
	vector<vector<unsigned int> > pressureIndex;
	vector<vector<vector<doubflo> > > velocityQs;
	vector<vector<unsigned int> > velocityIndex;
	vector<vector<vector<doubflo> > > outflowQs;
	vector<vector<unsigned int> > outflowIndex;
	//Geom-Daten werden mit dem Methode initArray() auf die Arrays kopiert. kein GetSet!

	//cout << "5: MyID: " << para->getMyID() << endl;
	noslipQs.resize(level+1);
	slipQs.resize(level+1);
	pressureQs.resize(level+1);
	velocityQs.resize(level+1);
	outflowQs.resize(level+1);

	noslipIndex.resize(level+1);
	slipIndex.resize(level+1);
	pressureIndex.resize(level+1);
	velocityIndex.resize(level+1);
	outflowIndex.resize(level+1);

	//cout << "6: MyID: " << para->getMyID() << endl;

	for(int i=0; i<=level;i++) {
		noslipQs[i].resize(27);
		slipQs[i].resize(27);
		pressureQs[i].resize(27);
		velocityQs[i].resize(27);
		outflowQs[i].resize(27);
	}

	//cout << "7: MyID: " << para->getMyID() << endl;

	for (int i = 0; i < t; i++){
		if (this->way[i] == "noSlip") { noslipQs = BC_Qs[i]->setgetBoundarys(noslipQs); noslipIndex = BC_Qs[i]->setgetIndex(noslipIndex); }
		else if (this->way[i] == "velocity") { velocityQs = BC_Qs[i]->setgetBoundarys(velocityQs); velocityIndex = BC_Qs[i]->setgetIndex(velocityIndex); }
		else if (this->way[i] == "pressure") { pressureQs = BC_Qs[i]->setgetBoundarys(pressureQs); pressureIndex = BC_Qs[i]->setgetIndex(pressureIndex); }
		else if (this->way[i] == "slip") { slipQs = BC_Qs[i]->setgetBoundarys(slipQs); slipIndex = BC_Qs[i]->setgetIndex(slipIndex); }
		else if (this->way[i] == "outflow") { outflowQs = BC_Qs[i]->setgetBoundarys(outflowQs); outflowIndex = BC_Qs[i]->setgetIndex(outflowIndex); }
	}



	//TEST Soeren---------------------------------------------------------------------------------------------------------------------------------------------------
	//if (TEST == true) {
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)velocityQs[i][0].size();
	//		cout << "Groesse VelocityQs lvl " << i << " : " << temp << endl;
	//	}
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)pressureQs[i][0].size();
	//		cout << "Groesse PressureQs lvl " << i << " : " << temp << endl;
	//	}
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)slipQs[i][0].size();
	//		cout << "Groesse slipQs lvl " << i << " : " << temp << endl;
	//	}
	//	for (int i = 0; i <= level; i++) {
	//		int temp = (int)noslipQs[i][0].size();
	//		cout << "Groesse noslipQs lvl " << i << " : " << temp << endl;
	//	}
	//	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//int ausgabe = 1;
	//	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	cout << endl;
	//	cout << "TEST PRESSURE Qs: " << endl;
	//	for (int j = 0; j <= level; j++) {
	//		if (pressureQs[j][0].size() > 0){
	//			for (int i = 0; i < ausgabe; i++) {
	//				cout << "PressureQs level " << j << ", eintrag " << i << " , index " << pressureIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << pressureQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//			for (unsigned int i = pressureQs[j][0].size() - ausgabe; i < pressureQs[j][0].size(); i++) {
	//				cout << "PressureQs level " << j << ", eintrag " << i << " , index " << pressureIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << pressureQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//		}
	//	}
	//	cout << endl;
	//	cout << "TEST VELOCITY Qs: " << endl;
	//	for (int j = 0; j <= level; j++) {
	//		if (velocityQs[j][0].size() > 0){
	//			for (int i = 0; i < ausgabe; i++) {
	//				cout << "velocityQs level " << j << ", eintrag " << i << " , index " << velocityIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << velocityQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//			for (unsigned int i = velocityQs[j][0].size() - ausgabe; i < velocityQs[j][0].size(); i++) {
	//				cout << "velocityQs level " << j << ", eintrag " << i << " , index " << velocityIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << velocityQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//		}
	//	}
	//	cout << endl;
		//cout << "TEST NOSLIP Qs: " << endl;
		//for (int j = 0; j <= level; j++) {
		//	if (noslipQs[j][0].size() > 0){
		//		for (int i = 0; i < ausgabe; i++) {
		//			cout << "noslipQs level " << j << ", eintrag " << i << " , index " << noslipIndex[j][i] << " :" << endl;
		//			for (int k = 0; k < 27; k++) {
		//				cout << noslipQs[j][k][i] << "|";
		//			}
		//			cout << endl;
		//			cout << endl;
		//		}
		//		for (unsigned int i = noslipQs[j][0].size() - ausgabe; i < noslipQs[j][0].size(); i++) {
		//			cout << "noslipQs level " << j << ", eintrag " << i << " , index " << noslipIndex[j][i] << " :" << endl;
		//			for (int k = 0; k < 27; k++) {
		//				cout << noslipQs[j][k][i] << "|";
		//			}
		//			cout << endl;
		//			cout << endl;
		//		}
		//	}
		//}
	//	cout << endl;
	//	cout << "TEST SLIP Qs: " << endl;
	//	for (int j = 0; j <= level; j++) {
	//		if (slipQs[j][0].size() > 0){
	//			for (int i = 0; i < ausgabe; i++) {
	//				cout << "slipQs level " << j << ", eintrag " << i << " , index " << slipIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << slipQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//			for (unsigned int i = slipQs[j][0].size() - ausgabe; i < slipQs[j][0].size(); i++) {
	//				cout << "slipQs level " << j << ", eintrag " << i << " , index " << slipIndex[j][i] << " :" << endl;
	//				for (int k = 0; k < 27; k++) {
	//					cout << slipQs[j][k][i] << "|";
	//				}
	//				cout << endl;
	//				cout << endl;
	//			}
	//		}
	//	}
	//}
	//ENDE TEST Soeren---------------------------------------------------------------------------------------------------------------------------------------------------
	
	//---------------------------------------------------------------------//AB HIER CODE FUER MARTIN////-------------------------------------------------------------------//
	//---------------------------------------------------------------------//
	//cout << "10: MyID: " << para->getMyID() << endl;

	for (int i = 0; i <= level; i++) {
		int temp1 = (int)pressureIndex[i].size();
		if (temp1 > 0)
		{
			cout << "Groesse Pressure:  " << i << " : " << temp1 << endl;
			//cout << "Groesse Pressure:  " << i << " : " << temp1 << "MyID: " << para->getMyID() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			doubflo* QQ = para->getParH(i)->QPress.q27[0];
			unsigned int sizeQ = para->getParH(i)->QPress.kQ;
			QforBoundaryConditions Q;
			Q.q27[dirE] = &QQ[dirE   *sizeQ];
			Q.q27[dirW] = &QQ[dirW   *sizeQ];
			Q.q27[dirN] = &QQ[dirN   *sizeQ];
			Q.q27[dirS] = &QQ[dirS   *sizeQ];
			Q.q27[dirT] = &QQ[dirT   *sizeQ];
			Q.q27[dirB] = &QQ[dirB   *sizeQ];
			Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
			Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
			Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
			Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
			Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
			Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
			Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
			Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
			Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
			Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
			Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
			Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
			Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
			Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
			Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
			Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
			Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
			Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
			Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
			Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
			Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = pressureQs.begin(); it != pressureQs.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = pressureIndex.begin(); it != pressureIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QPress.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// advection - diffusion stuff
			//cout << "vor advec diff" << endl;
			if (para->getDiffOn()==true){
				//////////////////////////////////////////////////////////////////////////
				//cout << "vor setzen von kTemp" << endl;
				para->getParH(i)->TempPress.kTemp = temp1;
				para->getParD(i)->TempPress.kTemp = temp1;
				cout << "Groesse TempPress.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
				//////////////////////////////////////////////////////////////////////////
				para->cudaAllocTempPressBC(i);
				//cout << "nach alloc" << endl;
				//////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp1; m++)
				{
					para->getParH(i)->TempPress.temp[m] = para->getTemperatureInit();
					para->getParH(i)->TempPress.velo[m] = (doubflo)0.0;
					para->getParH(i)->TempPress.k[m]    = para->getParH(i)->QPress.k[m];
				}
				//////////////////////////////////////////////////////////////////////////
				//cout << "vor copy" << endl;
				para->cudaCopyTempPressBCHD(i);
				//cout << "nach copy" << endl;
				//////////////////////////////////////////////////////////////////////////
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyPress(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife



	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp3 = (int)velocityIndex[i].size();
		if (temp3 > 0)
		{
			cout << "Groesse velocity level " << i << " : " << temp3 << endl;
			//cout << "Groesse velocity level:  " << i << " : " << temp3 << "MyID: " << para->getMyID() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			doubflo* QQ = para->getParH(i)->Qinflow.q27[0];
			unsigned int sizeQ = para->getParH(i)->Qinflow.kQ;
			QforBoundaryConditions Q;
			Q.q27[dirE] = &QQ[dirE   *sizeQ];
			Q.q27[dirW] = &QQ[dirW   *sizeQ];
			Q.q27[dirN] = &QQ[dirN   *sizeQ];
			Q.q27[dirS] = &QQ[dirS   *sizeQ];
			Q.q27[dirT] = &QQ[dirT   *sizeQ];
			Q.q27[dirB] = &QQ[dirB   *sizeQ];
			Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
			Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
			Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
			Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
			Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
			Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
			Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
			Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
			Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
			Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
			Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
			Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
			Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
			Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
			Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
			Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
			Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
			Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
			Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
			Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
			Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = velocityQs.begin(); it != velocityQs.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							Q.q27[j][n] = *it3;

							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = velocityIndex.begin(); it != velocityIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->Qinflow.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// advection - diffusion stuff
			if (para->getDiffOn()==true){
				//////////////////////////////////////////////////////////////////////////
				para->getParH(i)->TempVel.kTemp = temp3;
				para->getParD(i)->TempVel.kTemp = temp3;
				cout << "Groesse TempVel.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
				cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
				cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
				//////////////////////////////////////////////////////////////////////////
				para->cudaAllocTempVeloBC(i);
				//cout << "nach alloc " << endl;
				//////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp3; m++)
				{
					para->getParH(i)->TempVel.temp[m]      = para->getTemperatureInit();
					para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
					para->getParH(i)->TempVel.velo[m]      = para->getVelocity();
					para->getParH(i)->TempVel.k[m]         = para->getParH(i)->Qinflow.k[m];
				}
				//////////////////////////////////////////////////////////////////////////
				//cout << "vor copy " << endl;
				para->cudaCopyTempVeloBCHD(i);
				//cout << "nach copy " << endl;
				//////////////////////////////////////////////////////////////////////////
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyVeloBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getIsInflowNormal()){
				int temp = obj_inflowNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten InflowBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QInflowNormalX.kQ = temp;
					para->getParH(i)->QInflowNormalY.kQ = temp;
					para->getParH(i)->QInflowNormalZ.kQ = temp;
					para->getParD(i)->QInflowNormalX.kQ = para->getParH(i)->QInflowNormalX.kQ;
					para->getParD(i)->QInflowNormalY.kQ = para->getParH(i)->QInflowNormalY.kQ;
					para->getParD(i)->QInflowNormalZ.kQ = para->getParH(i)->QInflowNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocInflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_inflowNormalX->initIndex(para->getParH(i)->QInflowNormalX.k, i);
					obj_inflowNormalY->initIndex(para->getParH(i)->QInflowNormalY.k, i);
					obj_inflowNormalZ->initIndex(para->getParH(i)->QInflowNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					doubflo* QQX = para->getParH(i)->QInflowNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QInflowNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_inflowNormalX->initArray(QX.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Y
					doubflo* QQY = para->getParH(i)->QInflowNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QInflowNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_inflowNormalY->initArray(QY.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Z
					doubflo* QQZ = para->getParH(i)->QInflowNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QInflowNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_inflowNormalZ->initArray(QZ.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyInflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende oberste for schleife





	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp = (int)outflowIndex[i].size();
		if (temp > 0)
		{
			cout << "Groesse Outflow:  " << i << " : " << temp << endl;
			//cout << "Groesse Pressure:  " << i << " : " << temp1 << "MyID: " << para->getMyID() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			doubflo* QQ = para->getParH(i)->Qoutflow.q27[0];
			unsigned int sizeQ = para->getParH(i)->Qoutflow.kQ;
			QforBoundaryConditions Q;
			Q.q27[dirE] = &QQ[dirE   *sizeQ];
			Q.q27[dirW] = &QQ[dirW   *sizeQ];
			Q.q27[dirN] = &QQ[dirN   *sizeQ];
			Q.q27[dirS] = &QQ[dirS   *sizeQ];
			Q.q27[dirT] = &QQ[dirT   *sizeQ];
			Q.q27[dirB] = &QQ[dirB   *sizeQ];
			Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
			Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
			Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
			Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
			Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
			Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
			Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
			Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
			Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
			Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
			Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
			Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
			Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
			Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
			Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
			Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
			Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
			Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
			Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
			Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
			Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = outflowQs.begin(); it != outflowQs.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = outflowIndex.begin(); it != outflowIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->Qoutflow.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyOutflowBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getIsOutflowNormal()){
				int temp = obj_outflowNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten OutflowBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QOutflowNormalX.kQ = temp;
					para->getParH(i)->QOutflowNormalY.kQ = temp;
					para->getParH(i)->QOutflowNormalZ.kQ = temp;
					para->getParD(i)->QOutflowNormalX.kQ = para->getParH(i)->QOutflowNormalX.kQ;
					para->getParD(i)->QOutflowNormalY.kQ = para->getParH(i)->QOutflowNormalY.kQ;
					para->getParD(i)->QOutflowNormalZ.kQ = para->getParH(i)->QOutflowNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocOutflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_outflowNormalX->initIndex(para->getParH(i)->QOutflowNormalX.k, i);
					obj_outflowNormalY->initIndex(para->getParH(i)->QOutflowNormalY.k, i);
					obj_outflowNormalZ->initIndex(para->getParH(i)->QOutflowNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					doubflo* QQX = para->getParH(i)->QOutflowNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QOutflowNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_outflowNormalX->initArray(QX.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Y
					doubflo* QQY = para->getParH(i)->QOutflowNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QOutflowNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_outflowNormalY->initArray(QY.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Z
					doubflo* QQZ = para->getParH(i)->QOutflowNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QOutflowNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_outflowNormalZ->initArray(QZ.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyOutflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}//ende if
	}//ende oberste for schleife





	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp2 = (int)noslipQs[i][0].size();
		para->getParH(i)->QWall.kQ = temp2;
		para->getParD(i)->QWall.kQ = para->getParH(i)->QWall.kQ;
		para->getParH(i)->kQ = para->getParH(i)->QWall.kQ;
		para->getParD(i)->kQ = para->getParH(i)->QWall.kQ;
		if (temp2 > 0)
		{
			cout << "Groesse NoSlip: " << i << " : " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//para->getParH(i)->QWall.kQ = temp2;
			//para->getParD(i)->QWall.kQ = para->getParH(i)->QWall.kQ;
			//para->getParH(i)->kQ = para->getParH(i)->QWall.kQ;
			//para->getParD(i)->kQ = para->getParH(i)->QWall.kQ;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocWallBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			doubflo* QQ = para->getParH(i)->QWall.q27[0];
			unsigned int sizeQ = para->getParH(i)->QWall.kQ;
			QforBoundaryConditions Q;
			Q.q27[dirE] = &QQ[dirE   *sizeQ];
			Q.q27[dirW] = &QQ[dirW   *sizeQ];
			Q.q27[dirN] = &QQ[dirN   *sizeQ];
			Q.q27[dirS] = &QQ[dirS   *sizeQ];
			Q.q27[dirT] = &QQ[dirT   *sizeQ];
			Q.q27[dirB] = &QQ[dirB   *sizeQ];
			Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
			Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
			Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
			Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
			Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
			Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
			Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
			Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
			Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
			Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
			Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
			Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
			Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
			Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
			Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
			Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
			Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
			Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
			Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
			Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
			Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = noslipQs.begin(); it != noslipQs.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = noslipIndex.begin(); it != noslipIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QWall.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyWallBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}//ende oberste for schleife

	//--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp2 = (int)slipQs[i][0].size();
		if (temp2 > 0)
		{
			cout << "Groesse Slip: " << i << " : " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->getParH(i)->QSlip.kQ = temp2;
			para->getParD(i)->QSlip.kQ = para->getParH(i)->QSlip.kQ;
			para->getParH(i)->kSlipQ = para->getParH(i)->QSlip.kQ;
			para->getParD(i)->kSlipQ = para->getParH(i)->QSlip.kQ;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocSlipBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			doubflo* QQ = para->getParH(i)->QSlip.q27[0];
			unsigned int sizeQ = para->getParH(i)->QSlip.kQ;
			QforBoundaryConditions Q;
			Q.q27[dirE] = &QQ[dirE   *sizeQ];
			Q.q27[dirW] = &QQ[dirW   *sizeQ];
			Q.q27[dirN] = &QQ[dirN   *sizeQ];
			Q.q27[dirS] = &QQ[dirS   *sizeQ];
			Q.q27[dirT] = &QQ[dirT   *sizeQ];
			Q.q27[dirB] = &QQ[dirB   *sizeQ];
			Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
			Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
			Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
			Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
			Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
			Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
			Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
			Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
			Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
			Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
			Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
			Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
			Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
			Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
			Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
			Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
			Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
			Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
			Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
			Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
			Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<doubflo> > >::iterator it = slipQs.begin(); it != slipQs.end(); it++) {
				if (i == d) {
					for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
						for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = slipIndex.begin(); it != slipIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QSlip.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopySlipBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}//ende oberste for schleife



	//--------------------------------------------------------------------------//
	if (para->getIsGeo()){
		for (int i = 0; i <= level; i++) {
			int temp4 = obj_geomQ->getSize(i);
			para->getParH(i)->QGeom.kQ = temp4;
			para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
			if (temp4 > 0)
			{
				cout << "Groesse der Daten GeomBoundaryQs, Level " << i << " : " << temp4 << endl;
				cout << "Groesse der Daten GeomBoundaryQs, Level:  " << i << " : " << temp4 << "MyID: " << para->getMyID() << endl;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//para->getParH(i)->QGeom.kQ = temp4;
				//para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaAllocGeomBC(i);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				//////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_geomQ->initIndex(para->getParH(i)->QGeom.k, i);
				//////////////////////////////////////////////////////////////////////////
				//preprocessing
				doubflo* QQ = para->getParH(i)->QGeom.q27[0];
				unsigned int sizeQ = para->getParH(i)->QGeom.kQ;
				QforBoundaryConditions Q;
				Q.q27[dirE] = &QQ[dirE   *sizeQ];
				Q.q27[dirW] = &QQ[dirW   *sizeQ];
				Q.q27[dirN] = &QQ[dirN   *sizeQ];
				Q.q27[dirS] = &QQ[dirS   *sizeQ];
				Q.q27[dirT] = &QQ[dirT   *sizeQ];
				Q.q27[dirB] = &QQ[dirB   *sizeQ];
				Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
				Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
				Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
				Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
				Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
				Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
				Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
				Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
				Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
				Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
				Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
				Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
				Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
				Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
				Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
				Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
				Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
				Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
				Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
				Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
				Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
				//////////////////////////////////////////////////////////////////
				for (int j = 0; j<27; j++) {
					obj_geomQ->initArray(Q.q27[j], i, j);
				}//ende der for schleife
				//////////////////////////////////////////////////////////////////
				for(int test = 0; test < temp4; test++)
				{				
					Q.q27[dirZERO][test] = 0.0f;
				}
				//for(int test = 0; test < 3; test++)
				//{
				//	for (int tmp = 0; tmp < 27; tmp++)
				//	{
				//		cout <<"Kuhs: " << Q.q27[tmp][test]  << endl;				
				//	}
				//}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// advection - diffusion stuff
				if (para->getDiffOn()==true){
					//////////////////////////////////////////////////////////////////////////
					para->getParH(i)->Temp.kTemp = temp4;
					para->getParD(i)->Temp.kTemp = temp4;
					cout << "Groesse Temp.kTemp = " << para->getParH(i)->Temp.kTemp << endl;
					//////////////////////////////////////////////////////////////////////////
					para->cudaAllocTempNoSlipBC(i);
					//////////////////////////////////////////////////////////////////////////
					for (int m = 0; m < temp4; m++)
					{
						para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
						para->getParH(i)->Temp.k[m]    = para->getParH(i)->QGeom.k[m];
					}
					//////////////////////////////////////////////////////////////////////////
					para->cudaCopyTempNoSlipBCHD(i);
					//////////////////////////////////////////////////////////////////////////
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaCopyGeomBC(i);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			if (para->getIsGeoNormal()){
				int temp = obj_geomNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten GeomBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QGeomNormalX.kQ = temp;
					para->getParH(i)->QGeomNormalY.kQ = temp;
					para->getParH(i)->QGeomNormalZ.kQ = temp;
					para->getParD(i)->QGeomNormalX.kQ = para->getParH(i)->QGeomNormalX.kQ;
					para->getParD(i)->QGeomNormalY.kQ = para->getParH(i)->QGeomNormalY.kQ;
					para->getParD(i)->QGeomNormalZ.kQ = para->getParH(i)->QGeomNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocGeomNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_geomNormalX->initIndex(para->getParH(i)->QGeomNormalX.k, i);
					obj_geomNormalY->initIndex(para->getParH(i)->QGeomNormalY.k, i);
					obj_geomNormalZ->initIndex(para->getParH(i)->QGeomNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					doubflo* QQX = para->getParH(i)->QGeomNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QGeomNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_geomNormalX->initArray(QX.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Y
					doubflo* QQY = para->getParH(i)->QGeomNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QGeomNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_geomNormalY->initArray(QY.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing Z
					doubflo* QQZ = para->getParH(i)->QGeomNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QGeomNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j<27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_geomNormalZ->initArray(QZ.q27[j], i, j);
					}//ende der for schleife
					//////////////////////////////////////////////////////////////////

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyGeomNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}

			}
		}
	}//ende oberstes for
	//--------------------------------------------------------------------------//
	//---------------------------------------------------------------------//BIS HIER CODE FUER MARTIN////-------------------------------------------------------------------//


	for (int i = 0; i < t; i++) {
		delete BC_Qs[i];
	}

	delete BC_Qs;
	delete obj_geomQ;

	cout << "-----Ende BoundaryQs------" <<endl;
	
}

void Interface::setDimensions(Parameter* para)
{
	ifstream numberNodes;
	numberNodes.open(para->getnumberNodes().c_str(), ios::in);
	if (!numberNodes) {
		cerr << "Fehler beim Oeffnen von file_NumberNodes" << para->getnumberNodes() << endl;
		exit(1);
	}
	//Dimensionen einlesen
	string buffer;
	int bufferInt;
	std::vector<int> localGridNX;
	std::vector<int> localGridNY;
	std::vector<int> localGridNZ;

	for (/*unsigned*/ int i = 0; i <= para->getMaxLevel(); i++) {
		numberNodes >> buffer;
		numberNodes >> bufferInt;
		//cout << "X: " << bufferInt << endl;
		localGridNX.push_back(bufferInt);
		numberNodes >> bufferInt;
		//cout << "Y: " << bufferInt << endl;
		localGridNY.push_back(bufferInt);
		numberNodes >> bufferInt;
		//cout << "Z: " << bufferInt << endl;
		localGridNZ.push_back(bufferInt);
	}
	para->setGridX(localGridNX);
	para->setGridY(localGridNY);
	para->setGridZ(localGridNZ);
}

//Crap by Martin S.
void Interface::setBoundingBox(Parameter* para)
{
	ifstream numberNodes;
	numberNodes.open(para->getLBMvsSI().c_str(), ios::in);
	if (!numberNodes) {
		cerr << "Fehler beim Oeffnen von file_LBMvsSI" << endl;
		exit(1);
	}
	//Dimensionen Bounding Box einlesen
	doubflo bufferDoubflo;
	std::vector<doubflo> minX, maxX, minY, maxY, minZ, maxZ;

	for (int i = 0; i <= para->getMaxLevel(); i++) {
		numberNodes >> bufferDoubflo;
		//cout << "minX: " << bufferDoubflo << endl;
		minX.push_back(bufferDoubflo);
		numberNodes >> bufferDoubflo;
		//cout << "minY: " << bufferDoubflo << endl;
		minY.push_back(bufferDoubflo);
		numberNodes >> bufferDoubflo;
		//cout << "minZ: " << bufferDoubflo << endl;
		minZ.push_back(bufferDoubflo);
		numberNodes >> bufferDoubflo;
		//cout << "maxX: " << bufferDoubflo << endl;
		maxX.push_back(bufferDoubflo);
		numberNodes >> bufferDoubflo;
		//cout << "maxY: " << bufferDoubflo << endl;
		maxY.push_back(bufferDoubflo);
		numberNodes >> bufferDoubflo;
		//cout << "maxZ: " << bufferDoubflo << endl;
		maxZ.push_back(bufferDoubflo);
	}
	para->setMinCoordX(minX);
	para->setMinCoordY(minY);
	para->setMinCoordZ(minZ);
	para->setMaxCoordX(maxX);
	para->setMaxCoordY(maxY);
	para->setMaxCoordZ(maxZ);

	//much more crap by Martin S.
	std::string tmp;
	std::vector<doubflo> scaleLBMtoSI;
	std::vector<doubflo> translateLBMtoSI;

	while (getline(numberNodes, tmp))
	{
		if (tmp == "from LBM to SI")
		{
			numberNodes >> tmp;
			numberNodes >> bufferDoubflo;
			scaleLBMtoSI.push_back(bufferDoubflo);
			//cout << "Scale LBM to SI X: " << bufferDoubflo << endl;
			numberNodes >> bufferDoubflo;
			scaleLBMtoSI.push_back(bufferDoubflo);
			//cout << "Scale LBM to SI Y: " << bufferDoubflo << endl;
			numberNodes >> bufferDoubflo;
			scaleLBMtoSI.push_back(bufferDoubflo);
			//cout << "Scale LBM to SI Z: " << bufferDoubflo << endl;
			numberNodes >> tmp;
			numberNodes >> bufferDoubflo;
			translateLBMtoSI.push_back(bufferDoubflo);
			//cout << "Translate LBM to SI X: " << bufferDoubflo << endl;
			numberNodes >> bufferDoubflo;
			translateLBMtoSI.push_back(bufferDoubflo);
			//cout << "Translate LBM to SI Y: " << bufferDoubflo << endl;
			numberNodes >> bufferDoubflo;
			translateLBMtoSI.push_back(bufferDoubflo);
			//cout << "Translate LBM to SI Z: " << bufferDoubflo << endl;
		}
	}
	para->setScaleLBMtoSI(scaleLBMtoSI);
	para->setTranslateLBMtoSI(translateLBMtoSI);

	numberNodes.close();

}

//Funktion zum ersetzten der periodischen Nachbarn----------------------------------------------------------------------------------
void Interface::initPeriodicNeigh(vector<vector<vector<unsigned int> > > periodV, vector<vector<unsigned int> > periodIndex,  string way) {
	vector<unsigned int>neighVec;
	vector<unsigned int>indexVec;
	
	int zaehler = 0;

	for(unsigned int i=0; i<neighX->getLevel();i++) {
		if(way=="periodic_y"){
			neighVec = neighY->getVec(i);
			//cout << " Test 1 " << endl;
		} 
		else if(way=="periodic_x"){
			neighVec = neighX->getVec(i);
			//cout << " Test 2 " << endl;
		}
		else if(way=="periodic_z"){
			neighVec = neighZ->getVec(i);
			//cout << " Test 3 " << endl;
		}
		else {
			cout << "Falscher String in periodicValue" << endl;
			exit(1);
		}

		//cout << " Test 4 " << endl;

		//cout << " i: " << i << endl;
		////cout << "Laenge neighVec:" << neighVec.size() << endl;
		//cout << "Laenge periodIndex: " << periodIndex[i].size() << endl;
		//cout << "Laenge 1 periodV: " << periodV.size() << endl;
		//cout << "Laenge 2 periodV: " << periodV[i].size() << endl;
		//cout << "Laenge 3 periodV: " << periodV[i][0].size() << endl;

		for (vector<unsigned int>::iterator it = periodIndex[i].begin(); it != periodIndex[i].end(); it++) {
			if(periodV[i][0][zaehler] != 0) {
				//cout << "Hier schreibe ich im Nachbarn: " << neighVec[*it] << " an der Stelle " << *it << " das rein: " << periodV[i][0][zaehler] << endl;
				neighVec[*it]=periodV[i][0][zaehler];
				//cout << " zaehler 1: " << zaehler << endl;				
			}
			//cout << " zaehler 2: " << zaehler << endl;				

			zaehler++;
		}


		if(way=="periodic_y"){
			neighY->setVec(i, neighVec);
		} 
		else if(way=="periodic_x"){
			neighX->setVec(i, neighVec);
		}
		else if(way=="periodic_z"){
			neighZ->setVec(i, neighVec);
		}

	}
}


void Interface::sortSystem(BoundaryValues **BC_Values, Parameter *para, int t) {

	for (int i = 0; i < t; i++){
		if (system[i].compare("inlet") == 0){ BC_Values[i] = new BoundaryValues(para->getinletBcValues()); }
		if (system[i].compare("outlet") == 0){ BC_Values[i] = new BoundaryValues(para->getoutletBcValues()); }
		if (system[i].compare("back") == 0){ BC_Values[i] = new BoundaryValues(para->getbackBcValues()); }
		if (system[i].compare("front") == 0){ BC_Values[i] = new BoundaryValues(para->getfrontBcValues()); }
		if (system[i].compare("top") == 0){ BC_Values[i] = new BoundaryValues(para->gettopBcValues()); }
		if (system[i].compare("bottom") == 0){ BC_Values[i] = new BoundaryValues(para->getbottomBcValues());}
	}
}

void Interface::sortSystem(BoundaryQs **BC_Qs, Parameter *para, int t) {

	for (int i = 0; i < t; i++){
		if (system[i].compare("inlet") == 0){ BC_Qs[i] = new BoundaryQs(para->getinletBcQs(), false); }
		if (system[i].compare("outlet") == 0){ BC_Qs[i] = new BoundaryQs(para->getoutletBcQs(), false); }
		if (system[i].compare("back") == 0){ BC_Qs[i] = new BoundaryQs(para->getbackBcQs(), false); }
		if (system[i].compare("front") == 0){ BC_Qs[i] = new BoundaryQs(para->getfrontBcQs(), false); }
		if (system[i].compare("top") == 0){ BC_Qs[i] = new BoundaryQs(para->gettopBcQs(), false); }
		if (system[i].compare("bottom") == 0){ BC_Qs[i] = new BoundaryQs(para->getbottomBcQs(), false); }
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Interface::rearrangeGeometry(Parameter* para, int lev)
{
	//redefine fluid nodes
	for (int index = 0; index < para->getParH(lev)->size_Mat_SP; index++)
	{
		if (para->getParH(lev)->geoSP[index] == GEO_FLUID_OLD)
		{
			para->getParH(lev)->geoSP[index] = GEO_FLUID;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
