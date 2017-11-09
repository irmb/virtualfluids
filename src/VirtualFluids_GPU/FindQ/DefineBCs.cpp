#include "FindQ/DefineBCs.h"
#include "FindQ/FindQ.h"

void findPressQShip(Parameter* para)
{
	//x = begin (0)
	findKforQPressX0(para, para->getCoarse());
	para->cudaAllocPressX0(para->getCoarse());
	findQPressX0(para, para->getCoarse());
	para->cudaCopyPressX0(para->getCoarse());
	//x = end (1)
	findKforQPressX1(para, para->getCoarse());
	para->cudaAllocPressX1(para->getCoarse());
	findQPressX1(para, para->getCoarse());
	para->cudaCopyPressX1(para->getCoarse());
	//for (int lev = para->getFine(); lev >= para->getCoarse(); lev--)
	//{
	//	findKforQPressX1(para, lev);
	//	para->cudaAllocPressX1(lev);
	//	findQPressX1(para, lev);
	//	para->cudaCopyPressX1(lev);
	//}
}




void findQ27(Parameter* para)
{
   for (int lev = para->getFine(); lev >= para->getCoarse(); lev--)
   {
      findKforQ(para, lev);

      para->getParH(lev)->kQ       = para->getParH(lev)->QWall.kQ;
	  para->getParD(lev)->kQ       = para->getParH(lev)->QWall.kQ;
	  para->getParD(lev)->QWall.kQ = para->getParH(lev)->QWall.kQ;
      printf("kQ= %d\n", para->getParH(lev)->kQ);

	  para->cudaAllocWallBC(lev);

      findQ(para, lev);

	  para->getParH(lev)->kQ       = para->getParH(lev)->QWall.kQ;
	  para->getParD(lev)->kQ       = para->getParH(lev)->QWall.kQ;
	  para->getParD(lev)->QWall.kQ = para->getParH(lev)->QWall.kQ;
      printf("kQ= %d\n", para->getParH(lev)->kQ);

	  para->cudaCopyWallBC(lev);
   }
}




void findBC27(Parameter* para)
{                                      
   if ( para->getMyID() == 0)
   {
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Inflow
      findKforQInflow(para);

      para->getParH(para->getCoarse())->kInflowQ = para->getParH(para->getCoarse())->Qinflow.kQ;
	  para->getParD(para->getCoarse())->kInflowQ = para->getParH(para->getCoarse())->Qinflow.kQ;
      printf("kInflowQ= %d\n", para->getParH(para->getCoarse())->kInflowQ);

	  para->cudaAllocVeloBC(0); //level = 0

      findQInflow(para);

	  para->cudaCopyVeloBC(0); //level = 0
   }

   //...!!!...next if gets a block comment for a simple test... 
   //if (  para->getMyID() == para->getNumprocs() - 1)
   //{
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //Outflow
	  // findKforQOutflow(para);

	  // para->getParH(para->getCoarse())->kOutflowQ = para->getParH(para->getCoarse())->Qoutflow.kQ;
	  // para->getParD(para->getCoarse())->kOutflowQ = para->getParH(para->getCoarse())->Qoutflow.kQ;
	  // printf("kOutflowQ= %d\n", para->getParH(para->getCoarse())->kOutflowQ);

	  // para->cudaAllocPressBC();

	  // findQOutflow(para);

	  // para->cudaCopyPressBC();
   //}


   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   SOUND WAVE BCs at the solid nodes (walls)        //////////////////////////////////////////////////////////////////////////////////////////////
   //                      by                            //////////////////////////////////////////////////////////////////////////////////////////////
   //               Maddin Schlaffer                     //////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //findKforQSchlaff( para->getParH(para->getCoarse())->nx,      
   //                  para->getParH(para->getCoarse())->ny, 
   //                  para->getParH(para->getCoarse())->gridNX,  
   //                  para->getParH(para->getCoarse())->gridNY, 
   //                  para->getParH(para->getCoarse())->gridNZ,  
   //                  para->getParH(para->getCoarse())->geo,   
   //                  QnH,
   //                  QsH,
   //                  QeH,
   //                  QwH);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //// Allocate Host Memory
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////North
   //unsigned int mem_size_N_Q_k = sizeof(int)*QnH.kQ;
   //unsigned int mem_size_N_Q_q = sizeof(doubflo)*QnH.kQ;
   //kNQ = QnH.kQ;
   //printf("kNQ= %d\n",kNQ);
   //cudaHostMemoryAllocate((void**) &QnH.q27[0], para->getD3Qxx()*mem_size_N_Q_q );
   //cudaHostMemoryAllocate((void**) &QnH.k,                       mem_size_N_Q_k );
   //cudaHostMemoryAllocate((void**) &VxNH,                        mem_size_N_Q_q );
   //cudaHostMemoryAllocate((void**) &VyNH,                        mem_size_N_Q_q );
   //cudaHostMemoryAllocate((void**) &VzNH,                        mem_size_N_Q_q );
   //cudaHostMemoryAllocate((void**) &deltaVNH,                    mem_size_N_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////South
   //unsigned int mem_size_S_Q_k = sizeof(int)*QsH.kQ;
   //unsigned int mem_size_S_Q_q = sizeof(doubflo)*QsH.kQ;
   //kSQ = QsH.kQ;
   //printf("kSQ= %d\n",kSQ);
   //cudaHostMemoryAllocate((void**) &QsH.q27[0], para->getD3Qxx()*mem_size_S_Q_q );
   //cudaHostMemoryAllocate((void**) &QsH.k,                       mem_size_S_Q_k );
   //cudaHostMemoryAllocate((void**) &VxSH,                        mem_size_S_Q_q );
   //cudaHostMemoryAllocate((void**) &VySH,                        mem_size_S_Q_q );
   //cudaHostMemoryAllocate((void**) &VzSH,                        mem_size_S_Q_q );
   //cudaHostMemoryAllocate((void**) &deltaVSH,                    mem_size_S_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////East
   //unsigned int mem_size_E_Q_k = sizeof(int)*QeH.kQ;
   //unsigned int mem_size_E_Q_q = sizeof(doubflo)*QeH.kQ;
   //kEQ = QeH.kQ;
   //printf("kEQ= %d\n",kEQ);
   //cudaHostMemoryAllocate((void**) &QeH.q27[0], para->getD3Qxx()*mem_size_E_Q_q );
   //cudaHostMemoryAllocate((void**) &QeH.k,                       mem_size_E_Q_k );
   //cudaHostMemoryAllocate((void**) &VxEH,                        mem_size_E_Q_q );
   //cudaHostMemoryAllocate((void**) &VyEH,                        mem_size_E_Q_q );
   //cudaHostMemoryAllocate((void**) &VzEH,                        mem_size_E_Q_q );
   //cudaHostMemoryAllocate((void**) &deltaVEH,                    mem_size_E_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////West
   //unsigned int mem_size_W_Q_k = sizeof(int)*QwH.kQ;
   //unsigned int mem_size_W_Q_q = sizeof(doubflo)*QwH.kQ;
   //kWQ = QwH.kQ;
   //printf("kWQ= %d\n",kWQ);
   //cudaHostMemoryAllocate((void**) &QwH.q27[0], para->getD3Qxx()*mem_size_W_Q_q );
   //cudaHostMemoryAllocate((void**) &QwH.k,                       mem_size_W_Q_k );
   //cudaHostMemoryAllocate((void**) &VxWH,                        mem_size_W_Q_q );
   //cudaHostMemoryAllocate((void**) &VyWH,                        mem_size_W_Q_q );
   //cudaHostMemoryAllocate((void**) &VzWH,                        mem_size_W_Q_q );
   //cudaHostMemoryAllocate((void**) &deltaVWH,                    mem_size_W_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////find inflow Q's on coarse grid
   //findQSchlaff(  para->getParH(para->getCoarse())->nx,      
   //               para->getParH(para->getCoarse())->ny, 
   //               para->getParH(para->getCoarse())->gridNX,  
   //               para->getParH(para->getCoarse())->gridNY, 
   //               para->getParH(para->getCoarse())->gridNZ,  
   //               para->getParH(para->getCoarse())->geo,   
   //               para->getParH(para->getCoarse())->k, 
   //               kNQ, VxNH, VyNH, VzNH, deltaVNH, QnH.q27[0], QnH,
   //               kSQ, VxSH, VySH, VzSH, deltaVSH, QsH.q27[0], QsH,
   //               kEQ, VxEH, VyEH, VzEH, deltaVEH, QeH.q27[0], QeH,
   //               kWQ, VxWH, VyWH, VzWH, deltaVWH, QwH.q27[0], QwH);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //// Allocate Device Memory
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////North   
   //cudaMemoryAllocate((void**) &QnD.q27[0], para->getD3Qxx()* mem_size_N_Q_q );
   //cudaMemoryAllocate((void**) &QnD.k,                        mem_size_N_Q_k );
   //cudaMemoryAllocate((void**) &VxND,                         mem_size_N_Q_q );
   //cudaMemoryAllocate((void**) &VyND,                         mem_size_N_Q_q );
   //cudaMemoryAllocate((void**) &VzND,                         mem_size_N_Q_q );
   //cudaMemoryAllocate((void**) &deltaVND,                     mem_size_N_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////South
   //cudaMemoryAllocate((void**) &QsD.q27[0], para->getD3Qxx()* mem_size_S_Q_q );
   //cudaMemoryAllocate((void**) &QsD.k,                        mem_size_S_Q_k );
   //cudaMemoryAllocate((void**) &VxSD,                         mem_size_S_Q_q );
   //cudaMemoryAllocate((void**) &VySD,                         mem_size_S_Q_q );
   //cudaMemoryAllocate((void**) &VzSD,                         mem_size_S_Q_q );
   //cudaMemoryAllocate((void**) &deltaVSD,                     mem_size_S_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////East
   //cudaMemoryAllocate((void**) &QeD.q27[0], para->getD3Qxx()* mem_size_E_Q_q );
   //cudaMemoryAllocate((void**) &QeD.k,                        mem_size_E_Q_k );
   //cudaMemoryAllocate((void**) &VxED,                         mem_size_E_Q_q );
   //cudaMemoryAllocate((void**) &VyED,                         mem_size_E_Q_q );
   //cudaMemoryAllocate((void**) &VzED,                         mem_size_E_Q_q );
   //cudaMemoryAllocate((void**) &deltaVED,                     mem_size_E_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////West
   //cudaMemoryAllocate((void**) &QwD.q27[0], para->getD3Qxx()* mem_size_W_Q_q );
   //cudaMemoryAllocate((void**) &QwD.k,                        mem_size_W_Q_k );
   //cudaMemoryAllocate((void**) &VxWD,                         mem_size_W_Q_q );
   //cudaMemoryAllocate((void**) &VyWD,                         mem_size_W_Q_q );
   //cudaMemoryAllocate((void**) &VzWD,                         mem_size_W_Q_q );
   //cudaMemoryAllocate((void**) &deltaVWD,                     mem_size_W_Q_q );
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //// Copy Host Memory to Device
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////North   
   //cudaMemoryCopy(QnD.q27[0], QnH.q27[0], para->getD3Qxx()* mem_size_N_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(QnD.k,      QnH.k,                        mem_size_N_Q_k,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VxND,       VxNH,                         mem_size_N_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VyND,       VyNH,                         mem_size_N_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VzND,       VzNH,                         mem_size_N_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(deltaVND,   deltaVNH,                     mem_size_N_Q_q,  cudaMemcpyHostToDevice);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////South
   //cudaMemoryCopy(QsD.q27[0], QsH.q27[0], para->getD3Qxx()* mem_size_S_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(QsD.k,      QsH.k,                        mem_size_S_Q_k,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VxSD,       VxSH,                         mem_size_S_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VySD,       VySH,                         mem_size_S_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VzSD,       VzSH,                         mem_size_S_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(deltaVSD,   deltaVSH,                     mem_size_S_Q_q,  cudaMemcpyHostToDevice);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////East
   //cudaMemoryCopy(QeD.q27[0], QeH.q27[0], para->getD3Qxx()* mem_size_E_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(QeD.k,      QeH.k,                        mem_size_E_Q_k,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VxED,       VxEH,                         mem_size_E_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VyED,       VyEH,                         mem_size_E_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VzED,       VzEH,                         mem_size_E_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(deltaVED,   deltaVEH,                     mem_size_E_Q_q,  cudaMemcpyHostToDevice);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////West
   //cudaMemoryCopy(QwD.q27[0], QwH.q27[0], para->getD3Qxx()* mem_size_W_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(QwD.k,      QwH.k,                        mem_size_W_Q_k,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VxWD,       VxWH,                         mem_size_W_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VyWD,       VyWH,                         mem_size_W_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(VzWD,       VzWH,                         mem_size_W_Q_q,  cudaMemcpyHostToDevice);
   //cudaMemoryCopy(deltaVWD,   deltaVWH,                     mem_size_W_Q_q,  cudaMemcpyHostToDevice);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

