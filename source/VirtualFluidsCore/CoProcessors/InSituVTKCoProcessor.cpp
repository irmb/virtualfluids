#ifdef VF_VTK

#include "InSituVTKCoProcessor.h"
#include <LBMKernel.h>
#include <D3Q27ETBCProcessor.h>
#include <vector>
#include <string>

#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <iostream>
#include <fstream>

using namespace std;

InSituVTKCoProcessor::InSituVTKCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
InSituVTKCoProcessor::InSituVTKCoProcessor( SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& configFile, SPtr<LBMUnitConverter> conv ) : CoProcessor(grid, s), conv(conv)
{
   gridRank  = Communicator::getInstance()->getProcessID(); 
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   readConfigFile(configFile);

   blockVector.resize(maxInitLevel+1);

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);
   }

   //initialization of communicator
   contr = vtkSmartPointer<vtkSocketController>::New();
   contr->Initialize();

   comm = vtkSmartPointer<vtkSocketCommunicator>::New();

   // Establish connection
   if (!comm->ConnectTo(wHostname.c_str(), wPort))
   {
      cerr << "Client error: Could not connect to the server."<< endl;
      return;
   }
 
}
//////////////////////////////////////////////////////////////////////////
InSituVTKCoProcessor::~InSituVTKCoProcessor()
{
   comm->CloseConnection();
}
//////////////////////////////////////////////////////////////////////////
void InSituVTKCoProcessor::process( double step )
{
   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "InSituVTKCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void InSituVTKCoProcessor::collectData( double step )
{
   int istep = static_cast<int>(step);

   unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   points = vtkPoints::New();
   unstructuredGrid->SetPoints(points);
   arrays[0] = vtkSmartPointer<vtkDoubleArray>::New();
   arrays[0]->SetNumberOfComponents(1);
   arrays[0]->SetName( "Rho" );
   arrays[1] = vtkSmartPointer<vtkDoubleArray>::New();
   arrays[1]->SetNumberOfComponents(1);
   arrays[1]->SetName( "Vx" );
   arrays[2] = vtkSmartPointer<vtkDoubleArray>::New();
   arrays[2]->SetNumberOfComponents(1);
   arrays[2]->SetName( "Vy" );
   arrays[3] = vtkSmartPointer<vtkDoubleArray>::New();
   arrays[3]->SetNumberOfComponents(1);
   arrays[3]->SetName( "Vz" );
   arrays[4] = vtkSmartPointer<vtkDoubleArray>::New();
   arrays[4]->SetNumberOfComponents(1);
   arrays[4]->SetName( "Press" );

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      for(SPtr<Block3D> block : blockVector[level])
      {
         if (block)
         {
            addData(block);
         }
      }
   }

   unstructuredGrid->GetPointData()->AddArray(arrays[0]);
   unstructuredGrid->GetPointData()->AddArray(arrays[1]);
   unstructuredGrid->GetPointData()->AddArray(arrays[2]);
   unstructuredGrid->GetPointData()->AddArray(arrays[3]);
   unstructuredGrid->GetPointData()->AddArray(arrays[4]);
   unstructuredGrid->GetPointData()->SetScalars(arrays[1]);

   if (!comm->Send(&istep, 1, 1, 11))
   {
      cerr << "Client error: Error sending data." << endl;
      return;
   }

   if (!comm->Send(unstructuredGrid, 1, 9))
   {
      cerr << "Server error: Error sending data." << endl;
      return;
   }

   //vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   //writer->SetInput(unstructuredGrid);
   //writer->SetFileName("test.vtu");
   //writer->SetDataModeToAscii();
   //writer->Update();

   UBLOG(logINFO,"InSituVTKCoProcessor step: " << istep);
}
//////////////////////////////////////////////////////////////////////////
void InSituVTKCoProcessor::addData( SPtr<Block3D> block )
{
   UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
   double         dx           = grid->getDeltaX(block);

   SPtr<LBMKernel> kernel = block->getKernel();
   SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();          
   SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();     
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal vx1,vx2,vx3,rho;

   //knotennummerierung faengt immer bei 0 an!
   int SWB,SEB,NEB,NWB,SWT,SET,NET,NWT;

   //Funktionszeiger
   typedef void (*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/,LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);

   CalcMacrosFct calcMacros = NULL;

   if(block->getKernel()->getCompressible())
   {
      calcMacros = &D3Q27System::calcCompMacroscopicValues;
   }
   else
   {
      calcMacros = &D3Q27System::calcIncompMacroscopicValues;
   }

   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int maxX1 = (int)(distributions->getNX1());
   int maxX2 = (int)(distributions->getNX2());
   int maxX3 = (int)(distributions->getNX3());

   //int minX1 = 1;
   //int minX2 = 1;
   //int minX3 = 1;

   //int maxX1 = (int)(distributions->getNX1());
   //int maxX2 = (int)(distributions->getNX2());
   //int maxX3 = (int)(distributions->getNX3());

   //nummern vergeben und node vector erstellen + daten sammeln
   CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3,-1);
   maxX1 -= 2;
   maxX2 -= 2;
   maxX3 -= 2;

   SPtr<BoundaryConditions> bcPtr;
   int nr = points->GetNumberOfPoints();

   double x[3];

   for(size_t ix3=minX3; ix3<=maxX3; ix3++)
   {
      for(size_t ix2=minX2; ix2<=maxX2; ix2++)
      {
         for(size_t ix1=minX1; ix1<=maxX1; ix1++)
         {
            if(!bcArray->isUndefined(ix1,ix2,ix3) && !bcArray->isSolid(ix1,ix2,ix3))
            {
               int index = 0;

               x[0] = double(val<1>(org) - val<1>(nodeOffset) + ix1*dx);
               x[1] = double(val<2>(org) - val<2>(nodeOffset) + ix2*dx);
               x[2] = double(val<3>(org) - val<3>(nodeOffset) + ix3*dx);

               points->InsertPoint((vtkIdType)nr, x);

               nodeNumbers(ix1,ix2,ix3) = nr++;
               
               distributions->getDistribution(f, ix1, ix2, ix3);
               calcMacros(f,rho,vx1,vx2,vx3);
               double press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);

               if (UbMath::isNaN(rho) || UbMath::isInfinity(rho)) 
                  UB_THROW( UbException(UB_EXARGS,"rho is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
               //rho=999.0;
               if (UbMath::isNaN(press) || UbMath::isInfinity(press)) 
                  UB_THROW( UbException(UB_EXARGS,"press is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
               //press=999.0;
               if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1)) 
                  UB_THROW( UbException(UB_EXARGS,"vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
               //vx1=999.0;
               if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2)) 
                  UB_THROW( UbException(UB_EXARGS,"vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
               //vx2=999.0;
               if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3)) 
                  UB_THROW( UbException(UB_EXARGS,"vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
               //vx3=999.0;

               arrays[0]->InsertNextValue(rho * conv->getFactorDensityLbToW2());
               arrays[1]->InsertNextValue(vx1 * conv->getFactorVelocityLbToW2());
               arrays[2]->InsertNextValue(vx2 * conv->getFactorVelocityLbToW2());
               arrays[3]->InsertNextValue(vx3 * conv->getFactorVelocityLbToW2());
               arrays[4]->InsertNextValue(press * conv->getFactorPressureLbToW2());
            }
         }
      }
   }
   maxX1 -= 1;
   maxX2 -= 1;
   maxX3 -= 1;

   vtkIdType ptIds[8];
   //cell vector erstellen
   for(int ix3=minX3; ix3<=maxX3; ix3++)
   {
      for(int ix2=minX2; ix2<=maxX2; ix2++)
      {
         for(int ix1=minX1; ix1<=maxX1; ix1++)
         {
            if(   (ptIds[0]=SWB=nodeNumbers( ix1  , ix2,   ix3   )) >= 0
               && (ptIds[1]=SEB=nodeNumbers( ix1+1, ix2,   ix3   )) >= 0
               && (ptIds[2]=NWB=nodeNumbers( ix1  , ix2+1, ix3   )) >= 0
               && (ptIds[3]=NEB=nodeNumbers( ix1+1, ix2+1, ix3   )) >= 0
               && (ptIds[4]=SWT=nodeNumbers( ix1  , ix2,   ix3+1 )) >= 0
               && (ptIds[5]=SET=nodeNumbers( ix1+1, ix2,   ix3+1 )) >= 0
               && (ptIds[6]=NWT=nodeNumbers( ix1  , ix2+1, ix3+1 )) >= 0 
               && (ptIds[7]=NET=nodeNumbers( ix1+1, ix2+1, ix3+1 )) >= 0
               )
            {
               unstructuredGrid->InsertNextCell((int)VTK_VOXEL, (vtkIdType)8, ptIds);
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void InSituVTKCoProcessor::readConfigFile( const std::string& configFile )
{
   ifstream ifs;
   ifs.open (configFile, ifstream::in);
   if(!ifs) throw UbException(UB_EXARGS,"can not open "+configFile);
 
   string dummy;
   int wRank = 0;
   getline(ifs, dummy);
   int np = Communicator::getInstance()->getNumberOfProcesses();

   while (ifs.good())
   {
      getline(ifs, dummy, ';');
      getline(ifs, wIP, ';');
      getline(ifs, wHostname, ';');
      getline(ifs, dummy);
      wPort = stoi(dummy);
      if(wRank == gridRank) break;
      wRank++;
   }
   ifs.close();
}

//////////////////////////////////////////////////////////////////////////

#endif


