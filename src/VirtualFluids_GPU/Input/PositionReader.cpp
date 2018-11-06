#include "PositionReader.h"

#include "VirtualFluidsBasics/basics/utilities/UbFileInputASCII.h"

static const int E    = dirE;		   //static const int E    = 0;
static const int W    = dirW;		   //static const int W    = 1;
static const int N    = dirN;		   //static const int N    = 2;
static const int S    = dirS;		   //static const int S    = 3;
static const int T    = dirT;		   //static const int T    = 4;
static const int B    = dirB;		   //static const int B    = 5;
static const int NE   = dirNE;		   //static const int NE   = 6;
static const int SW   = dirSW;		   //static const int SW   = 7;
static const int SE   = dirSE;		   //static const int SE   = 8;
static const int NW   = dirNW;		   //static const int NW   = 9;
static const int TE   = dirTE;		   //static const int TE   = 10;
static const int BW   = dirBW;		   //static const int BW   = 11;
static const int BE   = dirBE;		   //static const int BE   = 12;
static const int TW   = dirTW;		   //static const int TW   = 13;
static const int TN   = dirTN;		   //static const int TN   = 14;
static const int BS   = dirBS;		   //static const int BS   = 15;
static const int BN   = dirBN;		   //static const int BN   = 16;
static const int TS   = dirTS;		   //static const int TS   = 17;
static const int TNE  = dirTNE;		   //static const int TNE  = 18;
static const int TNW  = dirTNW;		   //static const int TNW  = 19;
static const int TSE  = dirTSE;		   //static const int TSE  = 20;
static const int TSW  = dirTSW;		   //static const int TSW  = 21;
static const int BNE  = dirBNE;		   //static const int BNE  = 22;
static const int BNW  = dirBNW;		   //static const int BNW  = 23;
static const int BSE  = dirBSE;		   //static const int BSE  = 24;
static const int BSW  = dirBSW;		   //static const int BSW  = 25;
static const int ZERO = dirZERO;	   //static const int ZERO = 26;
								 
								 
static const int INV_E   = dirE;		   //= W;  
static const int INV_W   = dirW;		   //= E;  
static const int INV_N   = dirN;		   //= S;  
static const int INV_S   = dirS;		   //= N;  
static const int INV_T   = dirT;		   //= B;  
static const int INV_B   = dirB;		   //= T;  
static const int INV_NE  = dirNE;		   //= SW; 
static const int INV_SW  = dirSW;		   //= NE; 
static const int INV_SE  = dirSE;		   //= NW; 
static const int INV_NW  = dirNW;		   //= SE; 
static const int INV_TE  = dirTE;		   //= BW; 
static const int INV_BW  = dirBW;		   //= TE; 
static const int INV_BE  = dirBE;		   //= TW; 
static const int INV_TW  = dirTW;		   //= BE; 
static const int INV_TN  = dirTN;		   //= BS; 
static const int INV_BS  = dirBS;		   //= TN; 
static const int INV_BN  = dirBN;		   //= TS; 
static const int INV_TS  = dirTS;		   //= BN; 
static const int INV_TNE = dirTNE;		   //= BSW;
static const int INV_TNW = dirTNW;		   //= BSE;
static const int INV_TSE = dirTSE;		   //= BNW;
static const int INV_TSW = dirTSW;		   //= BNE;
static const int INV_BNE = dirBNE;		   //= TSW;
static const int INV_BNW = dirBNW;		   //= TSE;
static const int INV_BSE = dirBSE;		   //= TNW;
static const int INV_BSW = dirBSW;		   //= TNE;
static const int INV_ZERO= dirZERO;		   //= ZERO;

//static  const int INVDIR[ENDDIR+1];

const int INVDIR[] = {	INV_E,   
						INV_W,  
						INV_N,  
						INV_S,  
						INV_T,  
						INV_B,  
						INV_NE, 
						INV_SW, 
						INV_SE, 
						INV_NW,
						INV_TE, 
						INV_BW, 
						INV_BE, 
						INV_TW, 
						INV_TN, 
						INV_BS, 
						INV_BN, 
						INV_TS,				//alt             Bezug                 neu
						INV_BNE,			//INV_TNE = BSW  = 25 = dirTSW  = TSW = INV_BNE  
						INV_TSE,			//INV_TNW = BSE  = 24 = dirBNW  = BNW = INV_TSE
						INV_BSE,			//INV_TSE = BNW  = 23 = dirTNW  = TNW = INV_BSE	
						INV_TNW,			//INV_TSW = BNE  = 22 = dirBSE  = BSE = INV_TNW	
						INV_BNW,			//INV_BNE = TSW  = 21 = dirTSE  = TSE = INV_BNW
						INV_TSW,			//INV_BNW = TSE  = 20 = dirBNE  = BNE = INV_TSW	
						INV_BSW,			//INV_BSE = TNW  = 19 = dirTNE  = TNE = INV_BSW
						INV_ZERO,			//INV_BSW = TNE  = 18 = dirZERO = ZERO= INV_ZERO
						INV_TNE};			//INV_ZERO= ZERO = 26 = dirBSW  = BSW = INV_TNE

static const int        optionDigits = 2;  //--> 2 bits für secondary Option
static const long long  maxOptionVal = ( 1<<optionDigits ) - 1; //2^3-1 -> 7
float q[27]; 
long long noslipBoundaryFlags;		
																     

//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileForAlloc(const std::string fileName, Parameter* para)
{
   UbFileInputASCII in(fileName);
   int maxlevel = in.readInteger();
   in.readLine();

   for (int level = 0; level <= maxlevel; level++)
   {
	   //Host
	   para->getParH(level)->size_Mat_SP = in.readInteger();
	   para->getParH(level)->mem_size_real_SP    = sizeof(real     ) * para->getParH(level)->size_Mat_SP;
	   para->getParH(level)->mem_size_int_SP        = sizeof(unsigned int) * para->getParH(level)->size_Mat_SP;
	   //Device
	   para->getParD(level)->size_Mat_SP            = para->getParH(level)->size_Mat_SP;
	   para->getParD(level)->mem_size_int_SP        = sizeof(unsigned int) * para->getParD(level)->size_Mat_SP;
	   para->getParD(level)->mem_size_real_SP    = sizeof(real     ) * para->getParD(level)->size_Mat_SP;

	   in.readLine();
	   for(unsigned int u = /*1*/0; u<para->getParH(level)->size_Mat_SP; u++)
	   {
		   in.readInteger();
	   }
	   in.readLine();
   }
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFile(const std::string fileName, std::string Type, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	if (Type == "geoVec")
	{
		for (int level = 0; level <= maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u = 0; u<para->getParH(level)->size_Mat_SP; u++)
			{
				para->getParH(level)->geoSP[u] = in.readInteger();       
			}
			in.readLine();
		}
	} 
	else if (Type == "neighborX")
	{
		for (int level = 0; level <= maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u = 0; u<para->getParH(level)->size_Mat_SP; u++)
			{
				//para->getParH(level)->neighborZ_SP[u] = in.readInteger();
				para->getParH(level)->neighborX_SP[u] = in.readInteger();
			}
			in.readLine();
		}
	}
	else if (Type == "neighborY")
	{
		for (int level = 0; level <= maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u = 0; u<para->getParH(level)->size_Mat_SP; u++)
			{
				para->getParH(level)->neighborY_SP[u] = in.readInteger();
			}
			in.readLine();
		}
	}
	else if (Type == "neighborZ")
	{
		for (int level = 0; level <= maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u = /*1*/0; u<para->getParH(level)->size_Mat_SP; u++)
			{
				//para->getParH(level)->neighborX_SP[u] = in.readInteger();
				para->getParH(level)->neighborZ_SP[u] = in.readInteger();
			}
			in.readLine();
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileInterfaceForAlloc(const std::string fileName, std::string Type, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();
	unsigned int test = 0;

	for (int level = 0; level < maxlevel; level++)
	{
		if (Type == "CF")
		{
			para->getParH(level)->K_CF                = in.readInteger();
			para->getParD(level)->K_CF                = para->getParH(level)->K_CF;
			para->getParH(level)->intCF.kCF           = para->getParH(level)->K_CF;
			para->getParD(level)->intCF.kCF           = para->getParH(level)->K_CF;
			para->getParH(level)->mem_size_kCF        = sizeof(unsigned int) * para->getParH(level)->K_CF;
			para->getParD(level)->mem_size_kCF        = sizeof(unsigned int) * para->getParD(level)->K_CF;
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
			{
				test = in.readInteger();
				if (test>=para->getParH(level)->size_Mat_SP)
				{
					printf("CF: das geht doch nicht!!!");
				}
			}
		} 
		else if (Type == "FC")
		{
			para->getParH(level)->K_FC                = in.readInteger();
			para->getParD(level)->K_FC                = para->getParH(level)->K_FC;
			para->getParH(level)->intFC.kFC           = para->getParH(level)->K_FC;
			para->getParD(level)->intFC.kFC           = para->getParH(level)->K_FC;
			para->getParH(level)->mem_size_kFC        = sizeof(unsigned int) * para->getParH(level)->K_FC;
			para->getParD(level)->mem_size_kFC        = sizeof(unsigned int) * para->getParD(level)->K_FC;
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
			{
				test = in.readInteger();
				if (test>=para->getParH(level)->size_Mat_SP)
				{
					printf("FC: das geht doch nicht!!!");
				}
			}
		}

		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileInterface(const std::string fileName, std::string Type, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	if (Type == "CFC")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
			{
				para->getParH(level)->intCF.ICellCFC[u] = in.readInteger();       
			}
			in.readLine();
		}
	} 
	else if (Type == "CFF")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
			{
				para->getParH(level)->intCF.ICellCFF[u] = in.readInteger();
			}
			in.readLine();
		}
	}
	else if (Type == "FCC")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
			{
				para->getParH(level)->intFC.ICellFCC[u] = in.readInteger();
			}
			in.readLine();
		}
	}
	else if (Type == "FCF")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
			{
				para->getParH(level)->intFC.ICellFCF[u] = in.readInteger();
			}
			in.readLine();
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileInterfaceOffsetForAlloc(const std::string fileName, std::string Type, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		if (Type == "CF")
		{
			para->getParH(level)->K_CF                = in.readInteger();
			para->getParD(level)->K_CF                = para->getParH(level)->K_CF;
			para->getParH(level)->intCF.kCF           = para->getParH(level)->K_CF;
			para->getParD(level)->intCF.kCF           = para->getParH(level)->K_CF;
			para->getParH(level)->mem_size_kCF_off    = sizeof(real) * para->getParH(level)->K_CF;
			para->getParD(level)->mem_size_kCF_off    = sizeof(real) * para->getParD(level)->K_CF;
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
			{
				in.readDouble();
				in.readDouble();
				in.readDouble();
			}
		} 
		else if (Type == "FC")
		{
			para->getParH(level)->K_FC                = in.readInteger();
			para->getParD(level)->K_FC                = para->getParH(level)->K_FC;
			para->getParH(level)->intFC.kFC           = para->getParH(level)->K_FC;
			para->getParD(level)->intFC.kFC           = para->getParH(level)->K_FC;
			para->getParH(level)->mem_size_kFC_off    = sizeof(real) * para->getParH(level)->K_FC;
			para->getParD(level)->mem_size_kFC_off    = sizeof(real) * para->getParD(level)->K_FC;
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
			{
				//in.readInteger();
				in.readDouble();
				in.readDouble();
				in.readDouble();
			}
		}

		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileInterfaceOffset(const std::string fileName, std::string Type, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	if (Type == "CF")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
			{
				para->getParH(level)->offCF.xOffCF[u] = (real)in.readDouble();       
				para->getParH(level)->offCF.yOffCF[u] = (real)in.readDouble();       
				para->getParH(level)->offCF.zOffCF[u] = (real)in.readDouble();       
			}
			in.readLine();
		}
	} 
	else if (Type == "FC")
	{
		for (int level = 0; level < maxlevel; level++)
		{
			in.readInteger();
			in.readLine();
			for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
			{
				//para->getParH(level)->offFC.xOffFC[u] = 0.0;       
				//para->getParH(level)->offFC.yOffFC[u] = 0.0;       
				//para->getParH(level)->offFC.zOffFC[u] = 0.0;  
				//in.readInteger();
				para->getParH(level)->offFC.xOffFC[u] = (real)in.readDouble();       
				para->getParH(level)->offFC.yOffFC[u] = (real)in.readDouble();       
				para->getParH(level)->offFC.zOffFC[u] = (real)in.readDouble();       
			}
			in.readLine();
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileNoSlipBcForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->kQ                  = in.readInteger();
		para->getParD(level)->kQ                  = para->getParH(level)->kQ;
		para->getParH(level)->QWall.kQ            = para->getParH(level)->kQ;
		para->getParD(level)->QWall.kQ            = para->getParH(level)->kQ;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kQ; u++)
		{
			in.readInteger();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileNoSlipBcQreadForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->kQread                  = in.readInteger();
		para->getParD(level)->kQread                  = para->getParH(level)->kQread;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kQread; u++)
		{
			in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileNoSlipBcPos(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kQ; u++)
		{
			para->getParH(level)->QWall.k[u] = in.readInteger();       
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileNoSlipBcQs(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	real test;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kQread; u++)
		{
			test = (real)in.readFloat();
			//if (test <= 0. && test >= -1.)
			//{
			//	test = -0.5;
			//} 
			//else if (test >= 0. && test <= 1.)
			//{
			//	test = 0.5;
			//}
			//else
			//{
			//	test = 100.;
			//}
			
			////TÄST -> shit
			//para->getParH(level)->QWall.qread[u] = 0.5f;
			//orig
			para->getParH(level)->QWall.qread[u] = test;
			//para->getParH(level)->QWall.qread[u] = (real)in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileNoSlipBcValue(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kQ; u++)
		{
			para->getParH(level)->QWall.valueQ[u] = in.readLongLong();     
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::findQs(Parameter* para)
{
	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		//////////////////////////////////////////////////////////////////
		//preprocessing
		real* QQ                  = para->getParH(lev)->QWall.q27[0]; 
		unsigned int sizeQ           = para->getParH(lev)->kQ; 
		QforBoundaryConditions Q;
		Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
		Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
		Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
		Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
		Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
		Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
		Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
		Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
		Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
		Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
		Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
		Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
		Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
		Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
		Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
		Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
		Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
		Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
		Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
		Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
		Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
		Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
		Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
		Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
		Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
		Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
		Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
		//////////////////////////////////////////////////////////////////
		//set all to -1.0
		for(unsigned int u=0; u<para->getParH(lev)->kQ; u++)
		{
			for (int dir = dirE; dir<=dirBSW; dir++)
			{
				Q.q27[dir][u] = -1.0f;
			}
		}
		//////////////////////////////////////////////////////////////////
		//find q-directions
		unsigned int v = 0;
		for(unsigned int u=0; u<para->getParH(lev)->kQ; u++)
		{
			noslipBoundaryFlags = para->getParH(lev)->QWall.valueQ[u];
			//noSlipBoundaryValueVec[level].push_back(bc->getNoSlipBoundary());
			for(int dir=dirE;dir<=dirBSW;dir++)
			{
				if(( ( noslipBoundaryFlags>>(optionDigits*dir) ) & maxOptionVal ) != 0 /*hasNoSlipBoundaryFlag(dir)*/)
				{
					Q.q27[INVDIR[dir]][u] = para->getParH(lev)->QWall.qread[v];
					//noSlipBoundaryQsVec[level].push_back(bc->getQ(LBMD3Q27System::INVDIR[dir]));
					v++;
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileSlipBcForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->kSlipQ              = in.readInteger();
		para->getParD(level)->kSlipQ              = para->getParH(level)->kSlipQ;
		para->getParH(level)->QSlip.kQ            = para->getParH(level)->kSlipQ;
		para->getParD(level)->QSlip.kQ            = para->getParH(level)->kSlipQ;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kSlipQ; u++)
		{
			in.readInteger();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileSlipBcQreadForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->kSlipQread                  = in.readInteger();
		para->getParD(level)->kSlipQread                  = para->getParH(level)->kSlipQread;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kSlipQread; u++)
		{
			in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileSlipBcPos(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kSlipQ; u++)
		{
			para->getParH(level)->QSlip.k[u] = in.readInteger();       
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileSlipBcQs(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kSlipQread; u++)
		{
			para->getParH(level)->QSlip.qread[u] = (real)in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFileSlipBcValue(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kSlipQ; u++)
		{
			para->getParH(level)->QSlip.valueQ[u] = in.readLongLong();     
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::findSlipQs(Parameter* para)
{
	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		//////////////////////////////////////////////////////////////////
		//preprocessing
		real* QQ                  = para->getParH(lev)->QSlip.q27[0]; 
		unsigned int sizeQ           = para->getParH(lev)->kSlipQ; 
		QforBoundaryConditions Q;
		Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
		Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
		Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
		Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
		Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
		Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
		Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
		Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
		Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
		Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
		Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
		Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
		Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
		Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
		Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
		Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
		Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
		Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
		Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
		Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
		Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
		Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
		Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
		Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
		Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
		Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
		Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
		//////////////////////////////////////////////////////////////////
		//set all to -1.0
		for(unsigned int u=0; u<para->getParH(lev)->kSlipQ; u++)
		{
			for (int dir = dirE; dir<=dirBSW; dir++)
			{
				Q.q27[dir][u] = -1.0f;
			}
		}
		//////////////////////////////////////////////////////////////////
		//find q-directions
		unsigned int v = 0;
		for(unsigned int u=0; u<para->getParH(lev)->kSlipQ; u++)
		{
			noslipBoundaryFlags = para->getParH(lev)->QSlip.valueQ[u];
			for(int dir=dirE;dir<=dirBSW;dir++)
			{
				if(( ( noslipBoundaryFlags>>(optionDigits*dir) ) & maxOptionVal ) != 0 )
				{
					Q.q27[INVDIR[dir]][u] = para->getParH(lev)->QSlip.qread[v];
					v++;
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePressBcForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	int test = 1;
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		test = in.readInteger();
		//if (test==0)
		//{
		//	printf("testlevel = %d\n", level);
		//	//continue;
		//}
		para->getParH(level)->kPressQ              = test;
		para->getParD(level)->kPressQ              = para->getParH(level)->kPressQ;
		para->getParH(level)->QPress.kQ            = para->getParH(level)->kPressQ;
		para->getParD(level)->QPress.kQ            = para->getParH(level)->kPressQ;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kPressQ; u++)
		{
			in.readInteger();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePressBcQreadForAlloc(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	int test = 1;
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		test = in.readInteger();
		//if (test==0)
		//{
		//	printf("testlevel1 = %d\n", level);
		//	//continue;
		//}
		para->getParH(level)->kPressQread                  = test;
		para->getParD(level)->kPressQread                  = para->getParH(level)->kPressQread;
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kPressQread; u++)
		{
			in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePressBcPos(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kPressQ; u++)
		{
			if (u==0)
			{
				para->getParH(level)->QPress.k[u] = in.readInteger();//0; 
				/*in.readInteger();*/
			}
			else
			{
				para->getParH(level)->QPress.k[u] = in.readInteger();
			}
			//setRhoBC
			para->getParH(level)->QPress.RhoBC[u] = (real)0.f;
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePressBcQs(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kPressQread; u++)
		{
			para->getParH(level)->QPress.qread[u] = (real)in.readFloat();
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePressBcValue(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		in.readInteger();
		in.readLine();
		for(unsigned int u=0; u<para->getParH(level)->kPressQ; u++)
		{
			para->getParH(level)->QPress.valueQ[u] = in.readLongLong();     
		}
		in.readLine();
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::findPressQs(Parameter* para)
{
	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		//////////////////////////////////////////////////////////////////
		//preprocessing
		real* QQ                  = para->getParH(lev)->QPress.q27[0]; 
		unsigned int sizeQ           = para->getParH(lev)->kPressQ; 
		QforBoundaryConditions Q;
		Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
		Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
		Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
		Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
		Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
		Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
		Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
		Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
		Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
		Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
		Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
		Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
		Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
		Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
		Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
		Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
		Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
		Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
		Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
		Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
		Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
		Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
		Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
		Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
		Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
		Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
		Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
		//////////////////////////////////////////////////////////////////
		//set all to -1.0
		for(unsigned int u=0; u<para->getParH(lev)->kPressQ; u++)
		{
			for (int dir = dirE; dir<=dirBSW; dir++)
			{
				Q.q27[dir][u] = -1.0f;
			}
		}
		//////////////////////////////////////////////////////////////////
		//find q-directions
		unsigned int v = 0;
		for(unsigned int u=0; u<para->getParH(lev)->kPressQ; u++)
		{
			noslipBoundaryFlags = para->getParH(lev)->QPress.valueQ[u];
			for(int dir=dirE;dir<=dirBSW;dir++)
			{
				if(( ( noslipBoundaryFlags>>(optionDigits*dir) ) & maxOptionVal ) != 0 )
				{
					Q.q27[INVDIR[dir]][u] = para->getParH(lev)->QPress.qread[v];
					v++;
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePropellerCylinderForAlloc(Parameter* para)
{
	UbFileInputASCII in(para->getpropellerCylinder());
	int test = 0, count = 0;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->QPropeller.kQ = in.readInteger();
		para->getParD(level)->QPropeller.kQ = para->getParH(level)->QPropeller.kQ;
		in.readLine();
		if (level == para->getFine())
		{
			for(int u=0; u<para->getParH(level)->QPropeller.kQ; u++)
			{
				test = in.readInteger();
				if (para->getParH(level)->geoSP[test] == GEO_FLUID)
				{
					count++;
				}
				////////////////////////////////////////////////////////////////////////
				//for(unsigned int ix3=0; ix3<para->getParH(level)->nz; ix3++)
				//{
				//	for(unsigned int ix2=0; ix2<para->getParH(level)->ny; ix2++)
				//	{
				//		for(unsigned int ix1=0; ix1<para->getParH(level)->nx; ix1++)
				//		{
				//			unsigned int m = para->getParH(level)->nx*(para->getParH(level)->ny*ix3 + ix2) + ix1;
				//			if (para->getParH(level)->k[m] == test)
				//			{
				//				if(para->getParH(level)->geo[m] == 1)
				//				{
				//					count++;									
				//				}
				//			}
				//		}
				//	}
				//}
				//count++;
				////////////////////////////////////////////////////////////////////////
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		else
		{
			for(int u=0; u<para->getParH(level)->QPropeller.kQ; u++)
			{
				in.readInteger();
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		para->getParH(level)->QPropeller.kQ = count;
		para->getParD(level)->QPropeller.kQ = para->getParH(level)->QPropeller.kQ;
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePropellerCylinder(Parameter* para)
{
	UbFileInputASCII in(para->getpropellerCylinder());
	int test = 0, count = 0;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		int allnodes = in.readInteger();
		in.readLine();
		if (level == para->getFine())
		{
			for(int u=0; u<allnodes; u++)
			{
				test = in.readInteger();
				////////////////////////////////////////////////////////////////////////
				if (para->getParH(level)->geoSP[test] == GEO_FLUID)
				{
					para->getParH(level)->QPropeller.k[count] = test; 
					para->getParH(level)->QPropeller.Vx[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.Vy[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.Vz[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.RhoBC[count] = 0.0f;									
					count++;
				}
				else
				{
					in.readDouble();
					in.readDouble();
					in.readDouble();
				}
				//para->getParH(level)->QPropeller.k[count] = test; 
				//para->getParH(level)->QPropeller.Vx[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vy[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vz[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vx[count]	  = 0.07f;
				//para->getParH(level)->QPropeller.Vy[count]	  = 0.0f;
				//para->getParH(level)->QPropeller.Vz[count]	  = 0.0f;
				in.readLine();
			}
		} 
		else
		{
			for(int u=0; u<allnodes; u++)
			{
				in.readInteger(); 
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		printf("allnodes = %d, count = %d\n", allnodes, count);
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::definePropellerQs(Parameter* para)
{
	//////////////////////////////////////////////////////////////////
	//preprocessing
	real* QQ                  = para->getParH(para->getFine())->QPropeller.q27[0]; 
	unsigned int sizeQ           = para->getParH(para->getFine())->QPropeller.kQ; 
	QforBoundaryConditions Q;
	Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
	Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
	Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
	Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
	Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
	Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
	Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
	Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
	Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
	Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
	Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
	Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
	Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
	Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
	Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
	Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
	Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
	Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
	Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
	Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
	Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
	Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
	Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
	Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
	Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
	Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
	Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
	//////////////////////////////////////////////////////////////////
	for(int u=0; u<para->getParH(para->getFine())->QPropeller.kQ; u++)
	{
		for (int dir = dirE; dir<=dirBSW; dir++)
		{
			if ((dir==dirE)  || 
				(dir==dirNE) || (dir==dirSE) || (dir==dirTE) || (dir==dirBE) ||
				(dir==dirTNE)|| (dir==dirBNE)|| (dir==dirTSE)|| (dir==dirBSE))
			{
				Q.q27[dir][u] = 1.0f;
			} 
			else
			{
				Q.q27[dir][u] = -1.0f;
			}
		}
	}
	//////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readMeasurePoints( Parameter* para ) 
{
	UbFileInputASCII in(para->getmeasurePoints());
	int numberOfAllNodes = in.readInteger();
	in.readLine();
	int tempLevel;
	MeasurePoints tempMP;
	//printf("done, init the values...\n");
	for (int u = 0; u < numberOfAllNodes; u++)
	{
		tempMP.name = in.readString(); 		
		//printf("done, read the name...\n");
		tempMP.k = in.readInteger();
		//printf("done, read k...\n");
		tempLevel = in.readInteger();
		//printf("done, read level...\n");
		in.readLine();
		//printf("done, read the values...\n");
		para->getParH(tempLevel)->MP.push_back(tempMP);
		//printf("done, put it into a vector...\n");
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////


