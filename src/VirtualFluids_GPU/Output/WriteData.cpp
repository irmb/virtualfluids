#include "Output/WriteData.h"
#include "Output/ClogWriter.h"
#include "Output/LogWriter.hpp"
#include "Output/VtkSGWriter.hpp"
#include "Output/InterfaceDebugWriter.hpp"
#include "Output/UnstructuredGridWriter.hpp"
#include "Output/kFullWriter.hpp"
#include "Output/PosWriter.hpp"
#include "Output/PosVecIntWriter.hpp"
#include "Output/interfaceWriter.hpp"
#include "Output/OffsetWriter.hpp"
#include "Output/MeasurePointWriter.hpp"

#include "Parameter/Parameter.h"


////////////////////////////////////////////////////////////////////////////////
void writeInit(SPtr<Parameter> para)
{
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//copy Data to host
		para->cudaCopyPrint(lev);
		if (para->getCalcMedian()) para->cudaCopyMedianPrint(lev);
		if (para->getUseWale())    para->cudaCopyTurbulentViscosityDH(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const unsigned int numberOfParts = para->getParH(lev)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
		std::vector<std::string> fname;
		std::vector<std::string> fname_med;
		////2nd and 3rd Moments
		//std::vector<std::string> fname2ndMoments;
		//std::vector<std::string> fname3rdMoments;
		//std::vector<std::string> fnameHigherMoments;
		for (unsigned int i = 1; i <= numberOfParts; i++)
		{
			fname.push_back(para->getFName()+"_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(para->getTInit())+".vtk");
			fname_med.push_back(para->getFName()+"_bin_median_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(para->getTInit())+".vtk");
			//fname2ndMoments.push_back(para->getFName()+"_2ndMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(para->getTInit())+".vtk");
			//fname3rdMoments.push_back(para->getFName()+"_3rdMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(para->getTInit())+".vtk");
			//fnameHigherMoments.push_back(para->getFName()+"_HigherMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(para->getTInit())+".vtk");
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::string ffname_test = para->getFName()+"_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ ".vtk";
		//std::string ffname_bin  = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ "PartI.vtk";
		//std::string ffname_bin2 = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ "PartII.vtk";
		//std::string ffname_bin_Points = para->getFName()+"_Points_"+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ ".vtk";
		//std::string ffname_bin_med = para->getFName()+"_bin_med_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ ".vtk";
		//std::string ffname_bin_test = para->getFName()+"_bin_test_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit())+ ".vtk";
		//std::string ffname_bin_eff = para->getFName()+"_bin_eff_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(para->getTInit())+".vtk";
		//std::string ffname_bin_Points_eff = para->getFName()+"_Points_"+"_bin_eff_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(para->getTInit())+".vtk";

		//std::string ffname_bin  = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(para->getTInit()) + ".vtk";

		std::string fname_qs  = para->getFName()+"_Qs_Lev_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + ".vtk";

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////2nd and 3rd Moments
		////std::string fname2ndMoments = para->getFName()+"_2ndMoments_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(para->getTInit())+".vtk";
		////std::string fname3rdMoments = para->getFName()+"_3rdMoments_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(para->getTInit())+".vtk";
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (para->getDiffOn()==true)
		{
			//printf("vor writeUnstrucuredGridLTConc");
			UnstrucuredGridWriter::writeUnstrucuredGridLTConc(para.get(), lev, fname);
			//printf("nach writeUnstrucuredGridLTConc");
			//VtkSGWriter::writeVTKsgThS(para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//                           para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//                           para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_test,                para->getParH(lev)->geoSP,    
			//                           para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP, 
			//                           para->getVelocityRatio(),    para->getDensityRatio(),   
			//                           para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//                           para->getParH(lev)->dx,      para->getParH(lev)->k,      para->getParH(lev)->Conc);
		} 
		else
		{
			if (para->getUseWale())
			{
				UnstrucuredGridWriter::writeUnstrucuredGridLTwithTurbulentViscosity(para.get(), lev, fname);
			}
			//else if (para->getSimulatePorousMedia())
			//{
			//	UnstrucuredGridWriter::writeUnstrucuredGridPM(para, lev, fname);
			//}
			else
			{
				UnstrucuredGridWriter::writeUnstrucuredGridLT(para.get(), lev, fname);
			}

			//Debug
            InterfaceDebugWriter::writeInterfaceLinesDebugCF(para.get());
            InterfaceDebugWriter::writeInterfaceLinesDebugFC(para.get());
            InterfaceDebugWriter::writeInterfaceLinesDebugCFCneighbor(para.get());
            InterfaceDebugWriter::writeInterfaceLinesDebugCFFneighbor(para.get());
            InterfaceDebugWriter::writeInterfaceLinesDebugFCCneighbor(para.get());
            InterfaceDebugWriter::writeInterfaceLinesDebugFCFneighbor(para.get());
			InterfaceDebugWriter::writeNeighborXLinesDebug(para.get());
			InterfaceDebugWriter::writeNeighborYLinesDebug(para.get());
			InterfaceDebugWriter::writeNeighborZLinesDebug(para.get());

			//if (para->getParH(lev)->QGeom.kQ > 0)
			//{
			//	UnstrucuredGridWriter::writeQs(para, lev, fname_qs);
			//}

			////2nd and 3rd Moments
			//if (para->getCalc2ndOrderMoments())  UnstrucuredGridWriter::writeUnstrucuredGridEff2ndMomentsLT(para, lev, fname2ndMoments);
			//if (para->getCalc3rdOrderMoments())  UnstrucuredGridWriter::writeUnstrucuredGridEff3rdMomentsLT(para, lev, fname3rdMoments);
			//if (para->getCalcHighOrderMoments()) UnstrucuredGridWriter::writeUnstrucuredGridEffHigherMomentsLT(para, lev, fnameHigherMoments);


			//UnstrucuredGridWriter::writeUnstrucuredGrid(para, lev, ffname_bin, ffname_bin_Points);
			//UnstrucuredGridWriter::writeUnstrucuredGridEff(para, lev, ffname_bin_eff, ffname_bin_Points_eff);

			//UnstrucuredGridWriter::writeUnstrucuredGridAsciiEff(para, lev, ffname_bin, ffname_bin_Points);

			//VtkSGWriter::writeVTKsgSPbinTEST(ffname_bin_test);

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//2ndMoments
			//UnstrucuredGridWriter::writeUnstrucuredGridEff2ndMoments(para, lev, fname2ndMoments);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//VtkSGWriter::writeVTKsgSPbinAS( para->getParH(lev)->nx,      para->getParH(lev)->ny,        para->getParH(lev)->nz,        STARTOFFX, STARTOFFY, STARTOFFZ, 
			//							 para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY,    para->getParH(lev)->gridNZ, 
			//							 para->getParH(lev)->startz,  para->getParH(lev)->endz,      ffname_bin,                    para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//							 para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,     para->getParH(lev)->vz_SP,     para->getParH(lev)->rho_SP,
			//							 para->getParH(lev)->press_SP,							     
			//							 para->getVelocityRatio(),    para->getDensityRatio(),       
			//							 para->getParH(lev)->distX,   para->getParH(lev)->distY,     para->getParH(lev)->distZ,
			//							 para->getParH(lev)->dx,      para->getParH(lev)->coordX_SP, para->getParH(lev)->coordY_SP, para->getParH(lev)->coordZ_SP);
			////median
			//VtkSGWriter::writeVTKmedSPbinAS(    para->getParH(lev)->nx,			para->getParH(lev)->ny,			para->getParH(lev)->nz,         STARTOFFX, STARTOFFY, STARTOFFZ, 
			//								 para->getParH(lev)->gridNX,		para->getParH(lev)->gridNY,		para->getParH(lev)->gridNZ, 
			//								 para->getParH(lev)->startz,		para->getParH(lev)->endz,		ffname_bin_med,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//								 para->getParH(lev)->vx_SP_Med,     para->getParH(lev)->vy_SP_Med,  para->getParH(lev)->vz_SP_Med,  para->getParH(lev)->rho_SP_Med,
			//								 para->getParH(lev)->press_SP_Med,  1,
			//								 para->getVelocityRatio(),          para->getDensityRatio(),   
			//								 para->getParH(lev)->distX,         para->getParH(lev)->distY,      para->getParH(lev)->distZ,
			//								 para->getParH(lev)->dx,            para->getParH(lev)->coordX_SP,  para->getParH(lev)->coordY_SP,  para->getParH(lev)->coordZ_SP);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//VtkSGWriter::writeVTKsgSPbin(para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//							 para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//							 para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_bin,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,
			//							 para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP,
			//							 para->getParH(lev)->press_SP,
			//							 para->getVelocityRatio(),    para->getDensityRatio(),   
			//							 para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//							 para->getParH(lev)->dx,      para->getParH(lev)->k);
			////median
			//VtkSGWriter::writeVTKmedSPbin(  para->getParH(lev)->nx,			para->getParH(lev)->ny,			para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//								para->getParH(lev)->gridNX,		para->getParH(lev)->gridNY,		para->getParH(lev)->gridNZ, 
			//								para->getParH(lev)->startz,		para->getParH(lev)->endz,		ffname_bin_med,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//								para->getParH(lev)->vx_SP_Med,  para->getParH(lev)->vy_SP_Med,  para->getParH(lev)->vz_SP_Med,  para->getParH(lev)->rho_SP_Med,
			//								para->getParH(lev)->press_SP_Med, 1,
			//								para->getVelocityRatio(),    para->getDensityRatio(),   
			//								para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//								para->getParH(lev)->dx,      para->getParH(lev)->k);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//VtkSGWriter::writeVTKsgSP( para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//                           para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//                           para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_test,                para->getParH(lev)->geoSP,    
			//                           para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP,
			//                           para->getParH(lev)->press_SP,
			//                           para->getVelocityRatio(),    para->getDensityRatio(),   
			//                           para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//                           para->getParH(lev)->dx,      para->getParH(lev)->k);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////














////////////////////////////////////////////////////////////////////////////////
void writeTimestep(Parameter* para, unsigned int t)
{
	////////////////////////////////////////////////////////////////////////////////
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (para->getUseWale())    para->cudaCopyTurbulentViscosityDH(lev);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const unsigned int numberOfParts = para->getParH(lev)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
		std::vector<std::string> fname;
		std::vector<std::string> fname_med;
		//2nd and 3rd Moments
		std::vector<std::string> fname2ndMoments;
		std::vector<std::string> fname3rdMoments;
		std::vector<std::string> fnameHigherMoments;
		for (unsigned int i = 1; i <= numberOfParts; i++)
		{
			fname.push_back(para->getFName()+"_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(t)+".vtk");
			fname_med.push_back(para->getFName()+"_bin_median_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(t)+".vtk");
			fname2ndMoments.push_back(para->getFName()+"_2ndMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(t)+".vtk");
			fname3rdMoments.push_back(para->getFName()+"_3rdMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(t)+".vtk");
			fnameHigherMoments.push_back(para->getFName()+"_HigherMoments_bin_lev_"+StringUtil::toString<int>(lev)+"_ID_"+StringUtil::toString<int>(para->getMyID())+"_Part_"+StringUtil::toString<int>(i)+"_t_"+StringUtil::toString<int>(t)+".vtk");
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::string ffname_test = para->getFName()+"_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		//std::string ffname_bin  = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(t)+ "PartI.vtk";
		//std::string ffname_bin2 = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(t)+ "PartII.vtk";
		//std::string ffname_bin_Points = para->getFName()+"_Points_"+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		//std::string ffname_bin_eff = para->getFName()+"_bin_eff_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		//std::string ffname_bin_Points_eff = para->getFName()+"_Points_"+"_bin_eff_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		//std::string ffname_bin_med = para->getFName()+"_bin_med_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";

		//std::string ffname_bin  = para->getFName()+"_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID()) + "_" +StringUtil::toString<int>(t) + ".vtk";

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////2nd and 3rd Moments
		//std::string fname2ndMoments = para->getFName()+"_2ndMoments_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		//std::string fname3rdMoments = para->getFName()+"_3rdMoments_bin_"+StringUtil::toString<int>(lev)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+".vtk";
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (para->getDiffOn()==true)
		{
			//printf("vor writeUnstrucuredGridLTConc");
			UnstrucuredGridWriter::writeUnstrucuredGridLTConc(para, lev, fname);
			//printf("nach writeUnstrucuredGridLTConc");
			//VtkSGWriter::writeVTKsgThS(para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//                           para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//                           para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_test,                para->getParH(lev)->geoSP,    
			//                           para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP, 
			//                           para->getVelocityRatio(),    para->getDensityRatio(),   
			//                           para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//                           para->getParH(lev)->dx,      para->getParH(lev)->k,      para->getParH(lev)->Conc);
		} 
		else
		{
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//UnstrucuredGridWriter::writeUnstrucuredGrid(para, lev, ffname_bin, ffname_bin_Points);
			//UnstrucuredGridWriter::writeUnstrucuredGridEff(para, lev, ffname_bin_eff, ffname_bin_Points_eff);
			//UnstrucuredGridWriter::writeUnstrucuredGridBig(para, lev, ffname_bin, ffname_bin2);


			if (para->getUseWale())
			{
				UnstrucuredGridWriter::writeUnstrucuredGridLTwithTurbulentViscosity(para, lev, fname);
			}
			//else if (para->getSimulatePorousMedia())
			//{
			//	UnstrucuredGridWriter::writeUnstrucuredGridPM(para, lev, fname);
			//}
			else
			{
				UnstrucuredGridWriter::writeUnstrucuredGridLT(para, lev, fname);
			}

			//2nd and 3rd Moments
			if (para->getCalc2ndOrderMoments())  UnstrucuredGridWriter::writeUnstrucuredGridEff2ndMomentsLT(para, lev, fname2ndMoments);
			if (para->getCalc3rdOrderMoments())  UnstrucuredGridWriter::writeUnstrucuredGridEff3rdMomentsLT(para, lev, fname3rdMoments);
			if (para->getCalcHighOrderMoments()) UnstrucuredGridWriter::writeUnstrucuredGridEffHigherMomentsLT(para, lev, fnameHigherMoments);



			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getCalcMedian() && ((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()) && ((t%(unsigned int)para->getclockCycleForMP())==0))
			{
				//UnstrucuredGridWriter::writeUnstrucuredGridEffMedian(para, lev, ffname_bin_med);
				UnstrucuredGridWriter::writeUnstrucuredGridMedianLT(para, lev, fname_med);
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//2ndMoments
			//UnstrucuredGridWriter::writeUnstrucuredGridEff2ndMoments(para, lev, fname2ndMoments);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//VtkSGWriter::writeVTKsgSPbinAS( para->getParH(lev)->nx,      para->getParH(lev)->ny,        para->getParH(lev)->nz,        STARTOFFX, STARTOFFY, STARTOFFZ, 
			//							 para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY,    para->getParH(lev)->gridNZ, 
			//							 para->getParH(lev)->startz,  para->getParH(lev)->endz,      ffname_bin,                    para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//							 para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,     para->getParH(lev)->vz_SP,     para->getParH(lev)->rho_SP,
			//							 para->getParH(lev)->press_SP,							     
			//							 para->getVelocityRatio(),    para->getDensityRatio(),       
			//							 para->getParH(lev)->distX,   para->getParH(lev)->distY,     para->getParH(lev)->distZ,
			//							 para->getParH(lev)->dx,      para->getParH(lev)->coordX_SP, para->getParH(lev)->coordY_SP, para->getParH(lev)->coordZ_SP);
			//median
			//if (((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()))
			//{
			//VtkSGWriter::writeVTKmedSPbinAS(    para->getParH(lev)->nx,			para->getParH(lev)->ny,			para->getParH(lev)->nz,         STARTOFFX, STARTOFFY, STARTOFFZ, 
			//								 para->getParH(lev)->gridNX,		para->getParH(lev)->gridNY,		para->getParH(lev)->gridNZ, 
			//								 para->getParH(lev)->startz,		para->getParH(lev)->endz,		ffname_bin_med,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//								 para->getParH(lev)->vx_SP_Med,     para->getParH(lev)->vy_SP_Med,  para->getParH(lev)->vz_SP_Med,  para->getParH(lev)->rho_SP_Med,
			//								 para->getParH(lev)->press_SP_Med,  tdiff,
			//								 para->getVelocityRatio(),          para->getDensityRatio(),   
			//								 para->getParH(lev)->distX,         para->getParH(lev)->distY,      para->getParH(lev)->distZ,
			//								 para->getParH(lev)->dx,            para->getParH(lev)->coordX_SP,  para->getParH(lev)->coordY_SP,  para->getParH(lev)->coordZ_SP);
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//VtkSGWriter::writeVTKsgSPbin(para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//					   para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//					   para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_bin,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//					   para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP,
			//					   para->getParH(lev)->press_SP,
			//					   para->getVelocityRatio(),    para->getDensityRatio(),   
			//					   para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//					   para->getParH(lev)->dx,      para->getParH(lev)->k);
			//if (((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()))
			//{
			// //Median
			// real tdiff = (real)t - (real)t_prev;
			// VtkSGWriter::writeVTKmedSPbin(para->getParH(lev)->nx,           para->getParH(lev)->ny,         para->getParH(lev)->nz,         STARTOFFX, STARTOFFY, STARTOFFZ, 
			//									    para->getParH(lev)->gridNX,       para->getParH(lev)->gridNY,     para->getParH(lev)->gridNZ, 
			//									    para->getParH(lev)->startz,       para->getParH(lev)->endz,       ffname_bin_med,                 para->getParH(lev)->geoSP,  para->getParH(lev)->geo,    
			//									    para->getParH(lev)->vx_SP_Med,    para->getParH(lev)->vy_SP_Med,  para->getParH(lev)->vz_SP_Med,  para->getParH(lev)->rho_SP_Med,
			//									    para->getParH(lev)->press_SP_Med, tdiff,
			//									    para->getVelocityRatio(),         para->getDensityRatio(),   
			//									    para->getParH(lev)->distX,        para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//									    para->getParH(lev)->dx,           para->getParH(lev)->k);
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//            //VtkSGWriter::writeVTKsgSP( para->getParH(lev)->nx,      para->getParH(lev)->ny,     para->getParH(lev)->nz, STARTOFFX, STARTOFFY, STARTOFFZ, 
			//            //                           para->getParH(lev)->gridNX,  para->getParH(lev)->gridNY, para->getParH(lev)->gridNZ, 
			//            //                           para->getParH(lev)->startz,  para->getParH(lev)->endz,   ffname_test,                para->getParH(lev)->geoSP,    
			//            //                           para->getParH(lev)->vx_SP,   para->getParH(lev)->vy_SP,  para->getParH(lev)->vz_SP,  para->getParH(lev)->rho_SP,
			//            //                           para->getParH(lev)->press_SP,
			//            //                           para->getVelocityRatio(),    para->getDensityRatio(),
			//            //                           para->getParH(lev)->distX,   para->getParH(lev)->distY,  para->getParH(lev)->distZ,
			//            //                           para->getParH(lev)->dx,      para->getParH(lev)->k);
			//}
			//////////////////////////////////////////////////////////////////////////
			//output << "\n Write MeasurePoints at t = " << t << " (level = " << lev <<")\n";
			////for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
			////{
			//for(int j = 0; j < (int)para->getParH(lev)->MP.size(); j++)
			//{
			// MeasurePointWriter::writeMeasurePoints(para, lev, j, (int)t);
			//}
		}                                                  
		//////////////////////////////////////////////////////////////////////////
	}

}
////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////
void writeParticle(Parameter* para, unsigned int t)
{
	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		//////////////////////////////////////////////////////////////////////////
		//set filename
		std::string fname = para->getFName()+StringUtil::toString<int>(lev)+StringUtil::toString<int>(para->getMyID())+"_t_"+StringUtil::toString<int>(t)+"_Particles.vtk";
		//////////////////////////////////////////////////////////////////////////
		//write particles
		UnstrucuredGridWriter::writeUnstrucuredParticles(para, lev, fname);
		//////////////////////////////////////////////////////////////////////////
	}
}
////////////////////////////////////////////////////////////////////////////////
