#include "Particles/Particles.h"
//#include "Output/UnstructuredGridWriter.hpp"
#include "Output/FileWriter.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <stdio.h>

#include <fstream>
#include <sstream>
#include <random>

//#include <math.h>
//#include "LB.h"

//using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void allocParticles(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		//////////////////////////////////////////////////////////////////////////
		//all level store the same number of time steps
		para->getParH(lev)->plp.numberOfTimestepsParticles = para->getTOut() * (unsigned int)pow(2,para->getParticleBasicLevel());
		para->getParD(lev)->plp.numberOfTimestepsParticles = para->getTOut() * (unsigned int)pow(2,para->getParticleBasicLevel());
		//////////////////////////////////////////////////////////////////////////
		//all level store the same number of Particles
		para->getParH(lev)->plp.numberOfParticles = para->getNumberOfParticles();
		para->getParD(lev)->plp.numberOfParticles = para->getParH(lev)->plp.numberOfParticles;
		//////////////////////////////////////////////////////////////////////////
		//set memory size
		para->getParH(lev)->plp.memSizeTimestep      = sizeof(unsigned int)*(para->getParH(lev)->plp.numberOfTimestepsParticles);
		para->getParH(lev)->plp.memSizeID            = sizeof(unsigned int)*(para->getParH(lev)->plp.numberOfParticles);
		para->getParH(lev)->plp.memSizereal       = sizeof(real)*(para->getParH(lev)->plp.numberOfParticles);
		para->getParH(lev)->plp.memSizerealAll	 = sizeof(real)*(para->getParH(lev)->plp.numberOfParticles * para->getParH(lev)->plp.numberOfTimestepsParticles);
		para->getParH(lev)->plp.memSizeBool          = sizeof(bool)*(para->getParH(lev)->plp.numberOfParticles);
		para->getParH(lev)->plp.memSizeBoolBC        = sizeof(bool)*(para->getParH(lev)->geometryBC.numberOfBCnodes);//depends on Geometry!!!
		para->getParD(lev)->plp.memSizeTimestep      = para->getParH(lev)->plp.memSizeTimestep;
		para->getParD(lev)->plp.memSizeID            = para->getParH(lev)->plp.memSizeID;        
		para->getParD(lev)->plp.memSizereal       = para->getParH(lev)->plp.memSizereal;
		para->getParD(lev)->plp.memSizerealAll	 = para->getParH(lev)->plp.memSizerealAll;
		para->getParD(lev)->plp.memSizeBool          = para->getParH(lev)->plp.memSizeBool;
		para->getParD(lev)->plp.memSizeBoolBC        = para->getParH(lev)->plp.memSizeBoolBC;
		//////////////////////////////////////////////////////////////////////////
		//alloc particles
		cudaMemoryManager->cudaAllocParticles(lev);
		//////////////////////////////////////////////////////////////////////////
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initParticles(Parameter* para)
{
	para->setOutputCount((unsigned int) 0);
	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		//////////////////////////////////////////////////////////////////////////
		//set the time steps
		for (unsigned int t = 0; t < para->getParH(lev)->plp.numberOfTimestepsParticles; t++)
		{
			para->getParH(lev)->plp.timestep[t] = t;
		}
		//////////////////////////////////////////////////////////////////////////
		//set bool stuck to the wall
		for (unsigned int p = 0; p < para->getParH(lev)->plp.numberOfParticles; p++)
		{
			para->getParH(lev)->plp.stuck[p] = false;
		}
		//////////////////////////////////////////////////////////////////////////
		//set bool "hot wall"
		for (uint h = 0; h < para->getParH(lev)->geometryBC.numberOfBCnodes; h++)
		{
			if (para->getParH(lev)->coordinateX[para->getParH(lev)->geometryBC.k[h]] < para->getStartXHotWall() || 
				para->getParH(lev)->coordinateX[para->getParH(lev)->geometryBC.k[h]] > para->getEndXHotWall())
			{
				para->getParH(lev)->plp.hot[h] = false;
			}
			else
			{
				para->getParH(lev)->plp.hot[h] = true;
			}
		}
		//////////////////////////////////////////////////////////////////////////
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.timestep, para->getParH(lev)->plp.timestep, para->getParH(lev)->plp.memSizeTimestep, cudaMemcpyHostToDevice));
		//////////////////////////////////////////////////////////////////////////
		//test output
		int maximumT = para->getParH(lev)->plp.numberOfParticles * para->getParH(lev)->plp.numberOfTimestepsParticles;
		for (int i = 0; i < maximumT; i++)
		{
			para->getParH(lev)->plp.coordXlocal[i]   = (real)0.0;
			para->getParH(lev)->plp.coordYlocal[i]   = (real)0.0; 
			para->getParH(lev)->plp.coordZlocal[i]   = (real)0.0; 
			para->getParH(lev)->plp.coordXabsolut[i] = (real)0.0; 
			para->getParH(lev)->plp.coordYabsolut[i] = (real)0.0; 
			para->getParH(lev)->plp.coordZabsolut[i] = (real)0.0; 
			para->getParH(lev)->plp.veloX[i]         = (real)0.0; 
			para->getParH(lev)->plp.veloY[i]         = (real)0.0; 
			para->getParH(lev)->plp.veloZ[i]         = (real)0.0; 
		}
		//////////////////////////////////////////////////////////////////////////
		// real centerX  =  1.0f;					//uebergabeparameter
		real centerY  = 10.5f;					//uebergabeparameter
		real centerZ  = 10.5f;					//uebergabeparameter
		real diameter = 15.0f;//21.0f;					//uebergabeparameter
		unsigned int numberOfParticleSizes = 41;	//uebergabeparameter
		//////////////////////////////////////////////////////////////////////////
		//random
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 eng(rd()); // seed the generator
		std::uniform_int_distribution<> distr(1, numberOfParticleSizes); // define the range
		std::uniform_real_distribution<> distr1(centerY-(diameter/2.), centerY+(diameter/2.));
		std::uniform_real_distribution<> distr2(centerZ-(diameter/2.), centerZ+(diameter/2.));
		std::vector< int > radius;
		std::vector< double > yCoordVec;
		std::vector< double > zCoordVec;
		double dist;
		//////////////////////////////////////////////////////////////////////////
		real dx = (real)(1.0/pow(2,lev));
		std::vector< int > tempID;
		//////////////////////////////////////////////////////////////////////////

		for (unsigned int i = 0; i < para->getParH(lev)->plp.numberOfParticles; i++)
		{
			radius.push_back(distr(eng));
			yCoordVec.push_back(distr1(eng));
			zCoordVec.push_back(distr2(eng));

			dist = pow((zCoordVec[i]-centerZ),2) + pow((yCoordVec[i]-centerY),2);
			if(dist > pow((diameter/2),2))
			{
				zCoordVec[i] = sqrt(pow((diameter/2),2) - pow((yCoordVec[i]-centerY),2)) + centerZ;
			}
			para->getParH(lev)->plp.coordXabsolut[i] = (real)2.0; 
			para->getParH(lev)->plp.coordYabsolut[i] = (real)yCoordVec[i]; 
			para->getParH(lev)->plp.coordZabsolut[i] = (real)zCoordVec[i]; 

			// find IDs
			for (unsigned int ii = 0; ii < para->getParH(lev)->numberOfNodes; ii++)
			{
				if ((para->getParH(lev)->coordinateX[ii] <= para->getParH(lev)->plp.coordXabsolut[i]) &&
					((para->getParH(lev)->plp.coordXabsolut[i] - para->getParH(lev)->coordinateX[ii]) <= dx))
				{
					tempID.push_back(ii);
				}
			}

			//std::cout << "tempID[i]: " << tempID[i] << endl;
			// local paricle IDs
			para->getParH(lev)->plp.ID[i] = i;
		}
		//std::cout << "size tempID: " << tempID.size() << endl;
		//////////////////////////////////////////////////////////////////////////
		for (unsigned int i = 0; i < para->getParH(lev)->plp.numberOfParticles; i++)
		{
			for (std::size_t ii = 0; ii < tempID.size(); ii++)
			{
				if ((para->getParH(lev)->coordinateY[tempID[ii]] <= para->getParH(lev)->plp.coordYabsolut[i]) &&
					((para->getParH(lev)->plp.coordYabsolut[i] - para->getParH(lev)->coordinateY[tempID[ii]]) <= dx) &&
					(para->getParH(lev)->coordinateZ[tempID[ii]] <= para->getParH(lev)->plp.coordZabsolut[i]) &&
					((para->getParH(lev)->plp.coordZabsolut[i] - para->getParH(lev)->coordinateZ[tempID[ii]]) <= dx))
				{
					para->getParH(lev)->plp.cellBaseID[i] = tempID[ii];
					para->getParH(lev)->plp.coordXlocal[i] = (para->getParH(lev)->plp.coordXabsolut[i] - para->getParH(lev)->coordinateX[tempID[ii]]);
					para->getParH(lev)->plp.coordYlocal[i] = (para->getParH(lev)->plp.coordYabsolut[i] - para->getParH(lev)->coordinateY[tempID[ii]]);
					para->getParH(lev)->plp.coordZlocal[i] = (para->getParH(lev)->plp.coordZabsolut[i] - para->getParH(lev)->coordinateZ[tempID[ii]]);
				}
			}
			//std::cout << "ID: " << i << ", CellBaseID: " << para->getParH(lev)->plp.cellBaseID[i] 
			//<< ", localX: " << para->getParH(lev)->plp.coordXlocal[i] 
			//<< ", localY: " << para->getParH(lev)->plp.coordYlocal[i] 
			//<< ", localZ: " << para->getParH(lev)->plp.coordZlocal[i] << endl;
		}
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		//calc 1
		size_t in_var_max = 50;			//declearation of max size of input data values
		std::vector<double>sum(in_var_max);	//std::vector sum is used for calculation
		std::vector<int>c(in_var_max);		//This vector c stores the no.of counts that each particle size classes has occurred
		std::vector <int> c_temp1(in_var_max);
		static int c_temp2[50];

		//This loop is used to assign the no.of counts of the various particle size classes which are randomly generated
		for(unsigned int i = 0; i < para->getParH(lev)->plp.numberOfParticles; i++)
		{
			c_temp1[radius[i]-1]+=1;
		}//end of for loop

		int numberOfValues;				//variable to store the number of actual values given in a file
		std::vector<double>d_a(in_var_max);	//vector to store the diameter of the actual data given
		std::vector<double>q_a(in_var_max);	//vector q_a to store vol % of the actual data given
		std::vector<double>q(in_var_max);	//vector q to store vol % of the calculated data
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		//read file
		std::ifstream fin;					//object fin to read the file data
		fin.open("Actual_data.dat");	//opens the desired file name

		fin>>numberOfValues;

		//The following loop reads the data from the file as long as the file reaches the END OF LINE
		int j=0;
		while(!fin.eof() && j < numberOfValues)
		{
			fin>>(d_a[j]);
			fin>>(q_a[j]);
			j++;
		}
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		//calc 2
		//this loop calculates the desired (q%) for each particle size classes from the following counts of respective classes
		for(int i = 0; i < numberOfValues; i++)
		{
			c_temp2[i]+=c_temp1[i];
			c[i] = c_temp2[i];
			sum[i] = (q_a[i]/numberOfValues) * para->getParH(lev)->plp.numberOfParticles;//calculates the total q% for each particle size classes
			q[i] = sum[i]/c[i];//desired (q%) value,obtained by dividing with the no.of counts for each particle size classes
		}
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		//post processing
		//For loop;loops over every element of q
		for(int k=0; k < numberOfValues; k++)
		{
			int count=0;				//to count the no.of occurance for reach the acceptable values
			while(q[k]>q_a[k])			//checks whether the obtained results are greater than actual results - if greater 					performs the following operations
			{
				q[k] = sum[k]/c[k];
				if(q[k]<q_a[k]) break;
				c[k]+=1;				//Here I added the no.of count of the particular particle size class by 1
				q[k] = sum[k]/c[k];		//calc the (q%) again
				count++;				//increase the count variable by 1
			}

			while(count>0)				//This loops runs only when count becomes greater than 0
			{
				if(k!=numberOfValues && c[k+1]!=0)	//if condition checks for the non zero term and also not the end of vector
					c[k+1]-=1;			//this reduce the no.of count of neighbouring particle size class by 1 to counteract 					the increase in before while loop
				count--;				//decreasing count variable by 1
			}
		}//end of for loop
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.hot,           para->getParH(lev)->plp.hot,           para->getParH(lev)->plp.memSizeBoolBC,     cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.stuck,         para->getParH(lev)->plp.stuck,         para->getParH(lev)->plp.memSizeBool,       cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.ID,            para->getParH(lev)->plp.ID,            para->getParH(lev)->plp.memSizeID,         cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.cellBaseID,    para->getParH(lev)->plp.cellBaseID,    para->getParH(lev)->plp.memSizeID,         cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordXlocal,   para->getParH(lev)->plp.coordXlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordYlocal,   para->getParH(lev)->plp.coordYlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordZlocal,   para->getParH(lev)->plp.coordZlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordXabsolut, para->getParH(lev)->plp.coordXabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordYabsolut, para->getParH(lev)->plp.coordYabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordZabsolut, para->getParH(lev)->plp.coordZabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloX,         para->getParH(lev)->plp.veloX,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloY,         para->getParH(lev)->plp.veloY,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloZ,         para->getParH(lev)->plp.veloZ,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
		////////////////////////////////////////////////////////////////////////////
		////put some particles in the flow domain
		//InitParticlesDevice(para->getParD(lev)->coordX_SP,
		//					para->getParD(lev)->coordY_SP,
		//					para->getParD(lev)->coordZ_SP,
		//					para->getParD(lev)->plp.coordXlocal,
		//					para->getParD(lev)->plp.coordYlocal,
		//					para->getParD(lev)->plp.coordZlocal,
		//					para->getParD(lev)->plp.coordXabsolut,
		//					para->getParD(lev)->plp.coordYabsolut,
		//					para->getParD(lev)->plp.coordZabsolut,
		//					para->getParD(lev)->plp.veloX,
		//					para->getParD(lev)->plp.veloY,
		//					para->getParD(lev)->plp.veloZ,
		//					para->getParD(lev)->plp.randomLocationInit,
		//					para->getParD(lev)->plp.ID,
		//					para->getParD(lev)->plp.cellBaseID,
		//					para->getParD(lev)->geoSP,
		//					para->getParD(lev)->neighborX_SP,
		//					para->getParD(lev)->neighborY_SP,
		//					para->getParD(lev)->neighborZ_SP,
		//					para->getParD(lev)->neighborWSB_SP,
		//					lev,
		//					para->getParD(lev)->plp.numberOfParticles,
		//					para->getParD(lev)->size_Mat_SP,
		//					para->getParD(lev)->numberofthreads);
		////////////////////////////////////////////////////////////////////////////
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void propagateParticles(Parameter* para, unsigned int t)
{
	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		unsigned int tPLP = t - (para->getOutputCount() * para->getParH(lev)->plp.numberOfTimestepsParticles);
		//////////////////////////////////////////////////////////////////////////
		//put some particles in the flow domain
		MoveParticlesDevice(para->getParD(lev)->coordinateX,
							para->getParD(lev)->coordinateY,
							para->getParD(lev)->coordinateZ,
							para->getParD(lev)->plp.coordXlocal,
							para->getParD(lev)->plp.coordYlocal,
							para->getParD(lev)->plp.coordZlocal,
							para->getParD(lev)->plp.coordXabsolut,
							para->getParD(lev)->plp.coordYabsolut,
							para->getParD(lev)->plp.coordZabsolut,
							para->getParD(lev)->plp.veloX,
							para->getParD(lev)->plp.veloY,
							para->getParD(lev)->plp.veloZ,
							para->getParD(lev)->distributions.f[0],    
							para->getParD(lev)->omega, 
							para->getParD(lev)->plp.ID,
							para->getParD(lev)->plp.cellBaseID,
							para->getParD(lev)->typeOfGridNode,
							para->getParD(lev)->neighborX,
							para->getParD(lev)->neighborY,
							para->getParD(lev)->neighborZ,
							para->getParD(lev)->neighborInverse,
							lev,
							tPLP,
							para->getParD(lev)->plp.numberOfTimestepsParticles,
							para->getParD(lev)->plp.numberOfParticles,
							para->getParD(lev)->numberOfNodes,
							para->getParD(lev)->numberofthreads,
							para->getParD(lev)->isEvenTimestep);

	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void copyAndPrintParticles(Parameter* para, CudaMemoryManager* cudaMemoryManager, unsigned int t, bool isInit)
{
	//cout << endl << " t " << t << endl;

	//////////////////////////////////////////////////////////////////////////
	//count the number of outputs
	if (!isInit)
	{
		para->setOutputCount(para->getOutputCount() + (unsigned int)1);
	}
	//////////////////////////////////////////////////////////////////////////

	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		//////////////////////////////////////////////////////////////////////////
		//copy particles from device to host
		cudaMemoryManager->cudaCopyParticles(lev);
		//////////////////////////////////////////////////////////////////////////
		////test output
		//int maximumT = para->getParH(lev)->plp.numberOfParticles * para->getParH(lev)->plp.numberOfTimestepsParticles;
		//cout << " maximumT : " << maximumT << endl; 
		//for (int i = 0; i < maximumT; i++)
		//{
		//	cout << " i      : " << i << endl; 
		//	cout << " Xlocal : " << para->getParH(lev)->plp.coordXlocal[i] << endl; 
		//	cout << " Ylocal : " << para->getParH(lev)->plp.coordYlocal[i] << endl; 
		//	cout << " Zlocal : " << para->getParH(lev)->plp.coordZlocal[i] << endl; 
		//	cout << " Xglobal: " << para->getParH(lev)->plp.coordXabsolut[i] << endl; 
		//	cout << " Yglobal: " << para->getParH(lev)->plp.coordYabsolut[i] << endl; 
		//	cout << " Zglobal: " << para->getParH(lev)->plp.coordZabsolut[i] << endl; 
		//	cout << " vX     : " << para->getParH(lev)->plp.veloX[i] << endl; 
		//	cout << " vY     : " << para->getParH(lev)->plp.veloY[i] << endl; 
		//	cout << " vZ     : " << para->getParH(lev)->plp.veloZ[i] << endl; 
		//}
	}

	//////////////////////////////////////////////////////////////////////////
	//write particles
    //FileWriter().writeParticle(para, t);
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	if (!isInit)
	{
		for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
		{
			//////////////////////////////////////////////////////////////////////////
			int maximumT = (para->getParH(lev)->plp.numberOfParticles * para->getParH(lev)->plp.numberOfTimestepsParticles);
			int maximumP = para->getParH(lev)->plp.numberOfParticles;
			int j = (maximumT-maximumP);
			for (int i = 0; i < maximumP; i++, j++)
			{
				para->getParH(lev)->plp.coordXlocal[i]   = para->getParH(lev)->plp.coordXlocal[j]  ;
				para->getParH(lev)->plp.coordYlocal[i]   = para->getParH(lev)->plp.coordYlocal[j]  ; 
				para->getParH(lev)->plp.coordZlocal[i]   = para->getParH(lev)->plp.coordZlocal[j]  ; 
				para->getParH(lev)->plp.coordXabsolut[i] = para->getParH(lev)->plp.coordXabsolut[j]; 
				para->getParH(lev)->plp.coordYabsolut[i] = para->getParH(lev)->plp.coordYabsolut[j]; 
				para->getParH(lev)->plp.coordZabsolut[i] = para->getParH(lev)->plp.coordZabsolut[j]; 
				para->getParH(lev)->plp.veloX[i]         = para->getParH(lev)->plp.veloX[j]        ; 
				para->getParH(lev)->plp.veloY[i]         = para->getParH(lev)->plp.veloY[j]        ; 
				para->getParH(lev)->plp.veloZ[i]         = para->getParH(lev)->plp.veloZ[j]        ; 
			}
			//para->getParH(lev)->plp.coordXlocal[0]   = para->getParH(lev)->plp.coordXlocal[maximumT]  ;
			//para->getParH(lev)->plp.coordYlocal[0]   = para->getParH(lev)->plp.coordYlocal[maximumT]  ; 
			//para->getParH(lev)->plp.coordZlocal[0]   = para->getParH(lev)->plp.coordZlocal[maximumT]  ; 
			//para->getParH(lev)->plp.coordXabsolut[0] = para->getParH(lev)->plp.coordXabsolut[maximumT]; 
			//para->getParH(lev)->plp.coordYabsolut[0] = para->getParH(lev)->plp.coordYabsolut[maximumT]; 
			//para->getParH(lev)->plp.coordZabsolut[0] = para->getParH(lev)->plp.coordZabsolut[maximumT]; 
			//para->getParH(lev)->plp.veloX[0]         = para->getParH(lev)->plp.veloX[maximumT]        ; 
			//para->getParH(lev)->plp.veloY[0]         = para->getParH(lev)->plp.veloY[maximumT]        ; 
			//para->getParH(lev)->plp.veloZ[0]         = para->getParH(lev)->plp.veloZ[maximumT]        ; 
			//////////////////////////////////////////////////////////////////////////
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordXlocal,   para->getParH(lev)->plp.coordXlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordYlocal,   para->getParH(lev)->plp.coordYlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordZlocal,   para->getParH(lev)->plp.coordZlocal,   para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordXabsolut, para->getParH(lev)->plp.coordXabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordYabsolut, para->getParH(lev)->plp.coordYabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.coordZabsolut, para->getParH(lev)->plp.coordZabsolut, para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloX,         para->getParH(lev)->plp.veloX,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloY,         para->getParH(lev)->plp.veloY,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			checkCudaErrors( cudaMemcpy(para->getParD(lev)->plp.veloZ,         para->getParH(lev)->plp.veloZ,         para->getParH(lev)->plp.memSizerealAll, cudaMemcpyHostToDevice));
			////////////////////////////////////////////////////////////////////////////
		}
	}
	//////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////












////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rearrangeGeometry(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
	for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
		int counter1 = 0;
		int counter2 = 0;
		//////////////////////////////////////////////////////////////////////////
		//redefine fluid nodes
		for (uint index = 0; index < para->getParH(lev)->numberOfNodes; index++)
		{
			if (para->getParH(lev)->typeOfGridNode[index] == GEO_FLUID_OLD)
			{
				para->getParH(lev)->typeOfGridNode[index] = GEO_FLUID;
				counter1++;
			}
			counter2++;
		}
		printf("number of changed fluid nodes: %d \n", counter1);
		printf("total number of nodes: %d \n", counter2);
		//////////////////////////////////////////////////////////////////////////
		//store the index information of the BC nodes in the geometry array 
		for (uint index = 0; index < para->getParH(lev)->geometryBC.numberOfBCnodes; index++)
		{
			para->getParH(lev)->typeOfGridNode[para->getParH(lev)->geometryBC.k[index]] = index + OFFSET_BCsInGeo;
		}
		//////////////////////////////////////////////////////////////////////////
		//copy geometry nodes to the device
		cudaMemoryManager->cudaCopySP(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





