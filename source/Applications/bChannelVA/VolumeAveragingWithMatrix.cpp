#include "Postprocessing.h"
#include <vector>
#include <vtkDataArray.h>
#include <vtkPointLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkIntArray.h>
#include "MemoryUtil.h"
#include <string>


using namespace std;

vector<double> MATRIX;
vector<double> VMATRIX;
vector<int>    FMATRIX;
vector<int>    GMATRIX;
vector<double> LMATRIX;
vector<double> ZDATAVX;
vector<double> ZDATAVXX;
vector<double> ZDATAVYY;
vector<double> ZDATAVZZ;
vector<double> ZDATAC;

int NX1, NX2, NX3, NX4;
double DELTAX;
double ORG[3];

enum GEO{Solid, Fluid};
enum FLAG{Notexist, Exist};
enum VALUE1{DV, Vx, Vxx, Vyy, Vzz};
enum VALUE2{Flag, Geo};

//IndexerX4X3X2X1:
//x1-reihen "liegen am stueck" im speicher
//optimaler schleifendurchlauf
//for(alle X4)
//  for(alle X3)
//    for(alle X2)
//      for(alle X1)
int index(int x1, int x2, int x3, int x4)
{
   return NX1*(NX2*(NX3*x4+ x3)+x2)+x1;
}
//////////////////////////////////////////////////////////////////////////
int index2(int x1, int x2, int x3)
{
   return NX2 * ( NX3 * x3 + x2) + x1;
}
//////////////////////////////////////////////////////////////////////////
void get_node_indexes( double x[3], int level, int ix[3] )
{
   ix[0]=cint((x[0]-0.25*level-ORG[0])/DELTAX);
   ix[1]=cint((x[1]-0.25*level-ORG[1])/DELTAX);
   ix[2]=cint((x[2]-0.25*level-ORG[2])/DELTAX);
}
//////////////////////////////////////////////////////////////////////////
void get_node_indexes( double x[3], int ix[3] )
{
   ix[0]=cint((x[0]-ORG[0])/DELTAX);
   ix[1]=cint((x[1]-ORG[1])/DELTAX);
   ix[2]=cint((x[2]-ORG[2])/DELTAX);
}
//////////////////////////////////////////////////////////////////////////
void get_nodes_coordinates(int ix[3], int level, double x[3])
{
   x[0] = ORG[0] + (double)ix[2]*DELTAX + 0.25*level;
   x[1] = ORG[1] + (double)ix[2]*DELTAX + 0.25*level;
   x[2] = ORG[2] + (double)ix[2]*DELTAX + 0.25*level;
}
//////////////////////////////////////////////////////////////////////////
void create_matrix(std::string dataName, std::string dataNameG, std::string dataNameRMS, int numberOfLevels)
{
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_create = vtkSmartPointer<vtkTimerLog>::New();

   print("read data set from "+dataName+": start");
   timer_read->StartTimer();
   vtkDataSet* dataSet(ReadDataSet(dataName.c_str()));
   print("read data set from "+dataName+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   print("read data set from "+dataNameG+": start");
   timer_read->StartTimer();
   vtkDataSet* dataSetGeo(ReadDataSet(dataNameG.c_str()));
   print("read data set from "+dataNameG+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   print("read data set from "+dataNameRMS+": start");
   timer_read->StartTimer();
   vtkDataSet* dataSetRMS(ReadDataSet(dataNameRMS.c_str()));
   print("read data set from "+dataNameRMS+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   int numberOfPoints = dataSetGeo->GetNumberOfPoints();

   double bounds[6];
   dataSetGeo->GetBounds(bounds);

   ORG[0] = bounds[0];
   ORG[1] = bounds[2];
   ORG[2] = bounds[4];

   NX1 = 5;
   NX2 = int((bounds[1] - bounds[0]) / DELTAX) + 1;
   NX3 = int((bounds[3] - bounds[2]) / DELTAX) + 1;
   NX4 = int((bounds[5] - bounds[4]) / DELTAX) + 1;

   print("NX1xNX2xNX3xNX4 = "+toString(NX1)+"x"+toString(NX2)+"x"+toString(NX3)+"x"+toString(NX4));

   MATRIX.resize(NX1*NX2*NX3*NX4, 0.0);
   VMATRIX.resize(NX1*NX2*NX3*NX4, 0.0);
   FMATRIX.resize(NX2*NX3*NX4, 0);
   GMATRIX.resize(NX2*NX3*NX4, 0);
   ZDATAVX.resize(NX4, 0.0);
   ZDATAVXX.resize(NX4, 0.0);
   ZDATAVYY.resize(NX4, 0.0);
   ZDATAVZZ.resize(NX4, 0.0);

   vtkSmartPointer<vtkDataArray> vxArray = dataSet->GetPointData()->GetArray("Mvx");
   vtkSmartPointer<vtkDataArray> vxxArray = dataSetRMS->GetPointData()->GetArray("Vxx");
   vtkSmartPointer<vtkDataArray> vyyArray = dataSetRMS->GetPointData()->GetArray("Vyy");
   vtkSmartPointer<vtkDataArray> vzzArray = dataSetRMS->GetPointData()->GetArray("Vzz");
   vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
   vtkSmartPointer<vtkDataArray> geoArray = dataSetGeo->GetPointData()->GetArray("Geometry");

   double x[3];
   int   ix[3];

   print("create matrix: start");
   timer_create->StartTimer();

   int max_level = numberOfLevels - 1;

   for (int i = 0; i < numberOfPoints; i++)
   {
      dataSetGeo->GetPoint(i, x);
      vtkIdType id = dataSet->FindPoint(x);
      if (!((int)id < 0))
      {
         double l = lArray->GetTuple1(i);
         get_node_indexes(x, l, ix);
         double dx = (double)(1<<(max_level-(int)l));
         double dV = dx*dx*dx;
         MATRIX[index(DV, ix[0], ix[1], ix[2])] = dV;
         MATRIX[index(Vx, ix[0], ix[1], ix[2])] = vxArray->GetTuple1(id);
         MATRIX[index(Vxx, ix[0], ix[1], ix[2])] = vxxArray->GetTuple1(id);
         MATRIX[index(Vyy, ix[0], ix[1], ix[2])] = vyyArray->GetTuple1(id);
         MATRIX[index(Vzz, ix[0], ix[1], ix[2])] = vzzArray->GetTuple1(id);
         FMATRIX[index2(ix[0], ix[1], ix[2])] = Exist;
         if (geoArray->GetTuple1(i) == 0)
         {
            GMATRIX[index2(ix[0], ix[1], ix[2])] = Fluid;
         } 
         else
         {
            GMATRIX[index2(ix[0], ix[1], ix[2])] = Solid;
         }
      }
   }

   dataSet->Delete();
   dataSetGeo->Delete();
   dataSetRMS->Delete();

   print("create matrix: end");
   timer_create->StopTimer();
   print("create matrix time: "+toString(timer_create->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void write_matrix_to_image_file(std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output+".vti";

   print("write data set to "+vtkfilename+": start");
   timer_write->StartTimer();

   vtkImageData* image = vtkImageData::New();
   double delta = DELTAX;
   int extent[6];
   extent[0] = 0;
   extent[2] = 0;
   extent[4] = 0;
   extent[1] = NX2-1;
   extent[3] = NX3-1;
   extent[5] = NX4-1;
   image->SetExtent(extent);
   image->SetOrigin(ORG[0], ORG[1], ORG[2]);
   image->SetSpacing(delta, delta, delta);

   vtkDoubleArray* mvxArray = vtkDoubleArray::New();
   mvxArray->SetNumberOfComponents(1);
   mvxArray->SetName( "Mvx" );

   vtkDoubleArray* mvxxArray = vtkDoubleArray::New();
   mvxxArray->SetNumberOfComponents(1);
   mvxxArray->SetName( "Vxx" );

   vtkDoubleArray* mvyyArray = vtkDoubleArray::New();
   mvyyArray->SetNumberOfComponents(1);
   mvyyArray->SetName( "Vyy" );

   vtkDoubleArray* mvzzArray = vtkDoubleArray::New();
   mvzzArray->SetNumberOfComponents(1);
   mvzzArray->SetName( "Vzz" );

   vtkDoubleArray* dxArray = vtkDoubleArray::New();
   dxArray->SetNumberOfComponents(1);
   dxArray->SetName( "dx" );

   vtkDoubleArray* fArray = vtkDoubleArray::New();
   fArray->SetNumberOfComponents(1);
   fArray->SetName( "flag" );

   vtkDoubleArray* geoArray = vtkDoubleArray::New();
   geoArray->SetNumberOfComponents(1);
   geoArray->SetName( "geo" );

   for(int x3 = 0; x3 < NX4; x3++)
      for(int x2 = 0; x2 < NX3; x2++)
         for(int x1 = 0; x1 < NX2; x1++)   
         {
            dxArray->InsertNextValue(MATRIX[index(DV, x1, x2, x3)]);
            mvxArray->InsertNextValue(VMATRIX[index(Vx, x1, x2, x3)]);
            mvxxArray->InsertNextValue(VMATRIX[index(Vxx, x1, x2, x3)]);
            mvyyArray->InsertNextValue(VMATRIX[index(Vyy, x1, x2, x3)]);
            mvzzArray->InsertNextValue(VMATRIX[index(Vzz, x1, x2, x3)]);
            geoArray->InsertNextValue(GMATRIX[index2(x1, x2, x3)]);
            fArray->InsertNextValue(FMATRIX[index2(x1, x2, x3)]);
         }

         image->GetPointData()->AddArray(mvxArray);
         image->GetPointData()->AddArray(mvxxArray);
         image->GetPointData()->AddArray(mvyyArray);
         image->GetPointData()->AddArray(mvzzArray);
         image->GetPointData()->AddArray(dxArray);
         image->GetPointData()->AddArray(fArray);
         image->GetPointData()->AddArray(geoArray);

         vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
         writer->SetInputData(image);
         writer->SetFileName(vtkfilename.c_str());
         writer->SetDataModeToAscii();
         writer->Update();

         image->Delete();
         dxArray->Delete();
         mvxArray->Delete();
         mvxxArray->Delete();
         mvyyArray->Delete();
         mvzzArray->Delete();
         geoArray->Delete();
         fArray->Delete();
         writer->Delete();

         print("write data set: end");
         timer_write->StopTimer();
         print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void xy_averaging()
{
   print("xy averaging: start");
   vtkSmartPointer<vtkTimerLog> timer_sa = vtkSmartPointer<vtkTimerLog>::New();
   timer_sa->StartTimer();
   for(int x3 = 0; x3 < NX4; x3++)
   {
      double mvx = 0;
      double mvxx = 0;
      double mvyy = 0;
      double mvzz = 0;
      double i = 0;
      double volume = 0;

      for(int x2 = 0; x2 < NX3; x2++)
      {
         for(int x1 = 0; x1 < NX2; x1++)   
         {
            if (FMATRIX[index2(x1, x2, x3)] == Exist)
            {
               mvx += MATRIX[index(Vx, x1, x2, x3)];
               mvxx += MATRIX[index(Vxx, x1, x2, x3)];
               mvyy += MATRIX[index(Vyy, x1, x2, x3)];
               mvzz += MATRIX[index(Vzz, x1, x2, x3)];
               //volume += MATRIX[index(DV, x1, x2, x3)];
               i++;
            }
         }
      }
      if(i>0)
      {
         ZDATAVX[x3] = mvx / i ;/// (volume/i);
         ZDATAVXX[x3] = mvxx / i;
         ZDATAVYY[x3] = mvyy / i;
         ZDATAVZZ[x3] = mvzz / i;
      }
   }
   timer_sa->StopTimer();
   print("xy averaging: end");
   print("xy averaging time: "+toString(timer_sa->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void correct_index(int& x1, int& x2, int& x3)
{
   if (x1 < 0)   x1 = NX2 + x1;
   if (x1 >= NX2) x1 = x1 - NX2;

   if (x2 < 0)   x2 = NX3 + x2;
   if (x2 >= NX3) x2 = x2 - NX3;

   if (x3 < 0)   x3 = 0;
   if (x3 >= NX4) x3 = NX4 - 1;
}
//////////////////////////////////////////////////////////////////////////
void correct_index2(int& x1, int& x2, int& x3)
{
   if (x1 < 0)   x1 = 0;
   if (x1 >= NX2) x1 = NX2-1;

   if (x2 < 0)   x2 = 0;
   if (x2 >= NX3) x2 = NX3-1;

   if (x3 < 0)   x3 = 0;
   if (x3 >= NX4) x3 = NX4 - 1;
}
//////////////////////////////////////////////////////////////////////////
void v_averaging(double L)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();

   print("volume averaging: start");
   timer_averaging->StartTimer();

   //int X_Min = 0;
   //int X_Max = NX2;
   //int Y_Min = 0;
   //int Y_Max = NX3;
   //int Z_Min = 0;
   //int Z_Max = NX4;

   int x_min, x_max,y_min,y_max,z_min,z_max;
   //double dx,dy,dz;

   int l = int(L/2.0/DELTAX);

   x_min=l;
   x_max=l;
   y_min=l;
   y_max=l;
   z_min=l;
   z_max=l;

   double x_minC=L/2.0;
   double x_maxC=L/2.0;
   double y_minC=L/2.0;
   double y_maxC=L/2.0;
   double z_minC=L/2.0;
   double z_maxC=L/2.0;

   double lx_mean=(double)(x_minC+x_maxC)/2.0;
   double ly_mean=(double)(y_minC+y_maxC)/2.0;
   double lz_mean=(double)(z_minC+z_maxC)/2.0;
   //dy=dz=dx=1.0; 
   // = dx*dx*dx;

   double L3 = L*L*L;

   double norm = (lx_mean*lx_mean)*(ly_mean*ly_mean)*(lz_mean*lz_mean);
   int i = 0;

   timer_inloop->StartTimer();
   int p = 1000000;

#pragma omp parallel num_threads(8)
   {
#pragma omp for 
      for(int x3 = 0; x3 < NX4; x3++)
         for(int x2 = 0; x2 < NX3; x2++)
            for(int x1 = 0; x1 < NX2; x1++)   
            {
               if (FMATRIX[index(Flag, x1, x2, x3)] == Exist)
               {
                  if (i%p == 0 && i != 0)
                  {
                     timer_inloop->StartTimer();
                     print("point id = "+toString(i));
                  }

                  double dx3 = 0;
                  int j=0;
                  for(int z=x3-z_min;z<x3+z_max;z++)
                     for(int y=x2-y_min;y<x2+y_max;y++)
                        for(int x=x1-x_min;x<x1+x_max;x++)
                        {
                           int xx = x;
                           int yy = y;
                           int zz = z;
                           correct_index(xx, yy,zz);
                           if (FMATRIX[index(Flag,xx, yy, zz)] == Exist) 
                           {
                              double dx = MATRIX[index(DV, xx, yy, zz)];
                              dx3 += dx*dx*dx;
                              j++;
                           }
                        }

                        double factor = L3/dx3; 

                        for(int z=x3-z_min;z<x3+z_max;z++)
                           for(int y=x2-y_min;y<x2+y_max;y++)
                              for(int x=x1-x_min;x<x1+x_max;x++)
                              {
                                 int xx = x;
                                 int yy = y;
                                 int zz = z;
                                 correct_index(xx, yy,zz);
                                 if (FMATRIX[index(Flag, xx, yy, zz)] == Exist)
                                 {
                                    double dx = MATRIX[index(DV, xx, yy, zz)];
                                    double ddx3 = dx*dx*dx;
                                    double xc=ORG[0] + (double)x*DELTAX;
                                    double yc=ORG[1] + (double)y*DELTAX;
                                    double zc=ORG[2] + (double)z*DELTAX;
                                    double x1c=ORG[0] + (double)x1*DELTAX;
                                    double x2c=ORG[1] + (double)x2*DELTAX;
                                    double x3c=ORG[2] + (double)x3*DELTAX;
                                    //double mdv = m((double)x+dx/2.0-(double)x1,(double)y+dx/2.0-(double)x2,(double)z+dx/2.0-(double)x3,
                                    //               (double)x_min,(double)x_max,(double)y_min,(double)y_max,(double)z_min,(double)z_max,norm)*ddx3/factor;
                                    double mdv = m(xc+dx/2.0-x1c,yc+dx/2.0-x2c,zc+dx/2.0-x3c,x_minC,x_maxC,y_minC,y_maxC,z_minC,z_maxC,norm)*ddx3*factor;
                                    //double mdv = m(x+0.5-x1,y+0.5-x2,z+0.5-x3,x_min,x_max,y_min,y_max,z_min,z_max,norm)*factor;
                                    double gamma = (double)FMATRIX[index(Geo, xx, yy, zz)];
                                    VMATRIX[index(Vx, x1, x2, x3)]  += mdv*MATRIX[index(Vx, xx, yy, zz)]*gamma;
                                    VMATRIX[index(Vxx, x1, x2, x3)] += mdv*MATRIX[index(Vxx, xx, yy, zz)]*gamma;
                                    VMATRIX[index(Vyy, x1, x2, x3)] += mdv*MATRIX[index(Vyy, xx, yy, zz)]*gamma;
                                    VMATRIX[index(Vzz, x1, x2, x3)] += mdv*MATRIX[index(Vzz, xx, yy, zz)]*gamma;
                                 }
                              }


                              if (i%p == 0 && i != 0)
                              {
                                 timer_inloop->StopTimer();
                                 print("volume L3 = "+toString(L3));
                                 print("volume dx3 = "+toString(dx3));
                                 print("factor = "+toString(factor));
                                 print("number of points = "+toString(j));
                                 print("time per "+toString(p)+" points: "+toString(timer_inloop->GetElapsedTime())+" s");
                                 print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe()/1e9)+" GByte");
                                 timer_inloop->StartTimer();
                              }
                              i++;
               }
            }
   }
   timer_averaging->StopTimer();
   print("volume averaging: end");
   print("volume averaging time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
double G2(double x, double l)
{
   if(fabs(x)<=l) 
      return l-fabs(x);
   else
      return 0.0;
}
//////////////////////////////////////////////////////////////////////////
double m2(double x3, double l)
{
   return((G2(x3,l))/(l*l));
}
//////////////////////////////////////////////////////////////////////////
double gamma(int x1, int x2, int x3)
{
   if (FMATRIX[index(Geo, x1, x2, x3)] == Fluid)
   {
      return 1.0;
   } 
   else
   {
      return 0.0;
   }
}
//////////////////////////////////////////////////////////////////////////
void v_averaging2(double L)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();

   print("volume averaging: start");
   timer_averaging->StartTimer();

   print("L = "+toString(L));

   int l = cint(L/DELTAX);

   print("l = "+toString(l));

   int i = 0;

   timer_inloop->StartTimer();
   int p = 10000;
   double mdv;
   //#pragma omp parallel num_threads(8)
   {
      //#pragma omp for 
      for(int x3 = 0; x3 < NX4; x3++)
         for(int x2 = 0; x2 < NX3; x2++)
            for(int x1 = 0; x1 < NX2; x1++)   
            {
               if (FMATRIX[index2(x1, x2, x3)] == Exist)
               {
                  if (i%p == 0 && i != 0)
                  {
                     //timer_inloop->StartTimer();
                     print("point id = "+toString(i));
                  }

                  int vidx = index(DV, x1, x2, x3);

                  for(int z=-l;z<=+l;z++)
                     for(int y=-l;y<=+l;y++)
                        for(int x=-l;x<=+l;x++)
                        {
                           int xx = x1+x;
                           int yy = x2+y;
                           int zz = x3+z;
                           correct_index(xx, yy,zz);
                           int fgidx = index2(xx, yy, zz);

                           if (FMATRIX[fgidx] == Exist)
                           {
                              int midx = index(DV, xx, yy, zz);
                              double dV = MATRIX[midx];
                              mdv = m2((double)z, (double)l)*dV;
                              double g = (double)GMATRIX[fgidx];
                              VMATRIX[vidx+1] += mdv*MATRIX[midx+1]*g;
                              VMATRIX[vidx+2] += mdv*MATRIX[midx+2]*g;
                              VMATRIX[vidx+3] += mdv*MATRIX[midx+3]*g;
                              VMATRIX[vidx+4] += mdv*MATRIX[midx+4]*g;
                           }
                        }


                        if (i%p == 0 && i != 0)
                        {
                           //timer_inloop->StopTimer();
                           //print("volume L^3 = "+toString(L3));
                           //print("volume dV = "+toString(dV));
                           //print("factor = "+toString(factor));
                           //print("integral(m(x)) = "+toString(VMATRIX[index(Vx, x1, x2, x3)]));
                           //print("number of points = "+toString(j));
                           //print("k = "+toString(k));
                           //print("time per "+toString(p)+" points: "+toString(timer_inloop->GetElapsedTime())+" s");
                           print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe()/1e9)+" GByte");
                           //timer_inloop->StartTimer();
                        }
                        i++;
               }
            }
   }
   timer_averaging->StopTimer();
   print("volume averaging: end");
   print("volume averaging time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void s_averaging()
{
   vtkSmartPointer<vtkTimerLog> timer_sa = vtkSmartPointer<vtkTimerLog>::New();
   timer_sa->StartTimer();
   for(int x3 = 0; x3 < NX4; x3++)
   {
      double mvx = 0;
      double mvxx = 0;
      double mvyy = 0;
      double mvzz = 0;
      double i = 0;
      for(int x2 = 0; x2 < NX3; x2++)
      {
         for(int x1 = 0; x1 < NX2; x1++)   
         {
            if (FMATRIX[index2(x1, x2, x3)] == Exist && GMATRIX[index2(x1, x2, x3)] == Fluid)
            {
               mvx += VMATRIX[index(Vx, x1, x2, x3)];
               mvxx += VMATRIX[index(Vxx, x1, x2, x3)];
               mvyy += VMATRIX[index(Vyy, x1, x2, x3)];
               mvzz += VMATRIX[index(Vzz, x1, x2, x3)];
               i++;
            }
         }
      }
      if(i>0)
      {
         ZDATAVX[x3] = mvx / i;
         ZDATAVXX[x3] = mvxx / i;
         ZDATAVYY[x3] = mvyy / i;
         ZDATAVZZ[x3] = mvzz / i;
      }
   }
   timer_sa->StopTimer();
   print("surface averaging time: "+toString(timer_sa->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void write_ZDATA(std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();
   std::string fname = output+".dat";
   print("write z data to "+fname+": start");
   timer_write->StartTimer();
   std::ofstream ostr;

   ostr.open(fname.c_str(), std::ios_base::out);
   if(!ostr)
   { 
      ostr.clear();
   }

   ostr << "z" << "\t";
   ostr << "<Vx>" << "\t";
   ostr << "<Vxx>" << "\t";
   ostr << "<Vyy>" << "\t";
   ostr << "<Vzz>" << std::endl;

   for (int i = 0; i < (int)ZDATAVX.size(); i++)
   {
      if (ZDATAVX[i] && ZDATAVXX[i] && ZDATAVYY[i] && ZDATAVZZ[i] > 0)
      {
         //ostr << ORG[2] + (double)i*DELTAX << "\t";
         ostr << ZDATAC[i] << "\t";
         ostr << ZDATAVX[i] << "\t";
         ostr << ZDATAVXX[i] << "\t";
         ostr << ZDATAVYY[i] << "\t";
         ostr << ZDATAVZZ[i] << std::endl;
      }
   }

   ostr.close();

   print("write data set: end");
   timer_write->StopTimer();
   print("write z data time: "+toString(timer_write->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////
void read_image_data(std::string dataName)
{
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_create = vtkSmartPointer<vtkTimerLog>::New();

   print("read data set from "+dataName+": start");
   timer_read->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));
   print("read data set from "+dataName+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   vtkSmartPointer<vtkImageData> image = vtkImageData::SafeDownCast(dataSet);

   image->GetOrigin(ORG);

   int extent[6];
   image->GetExtent(extent);

   NX1 = 5;
   NX2 = extent[1]+1;
   NX3 = extent[3]+1;
   NX4 = extent[5]+1;

   print("NX1xNX2xNX3xNX4 = "+toString(NX1)+"x"+toString(NX2)+"x"+toString(NX3)+"x"+toString(NX4));

   MATRIX.resize(NX1*NX2*NX3*NX4, 0.0);
   VMATRIX.resize(NX1*NX2*NX3*NX4, 0.0);
   FMATRIX.resize(NX2*NX3*NX4, 0);
   GMATRIX.resize(NX2*NX3*NX4, 0);
   ZDATAVX.resize(NX4, 0.0);
   ZDATAVXX.resize(NX4, 0.0);
   ZDATAVYY.resize(NX4, 0.0);
   ZDATAVZZ.resize(NX4, 0.0);

   vtkSmartPointer<vtkDataArray> mvxArray = dataSet->GetPointData()->GetArray("Mvx");
   vtkSmartPointer<vtkDataArray> vxxArray = dataSet->GetPointData()->GetArray("Vxx");
   vtkSmartPointer<vtkDataArray> vyyArray = dataSet->GetPointData()->GetArray("Vyy");
   vtkSmartPointer<vtkDataArray> vzzArray = dataSet->GetPointData()->GetArray("Vzz");
   vtkSmartPointer<vtkDataArray> fArray = dataSet->GetPointData()->GetArray("flag");
   vtkSmartPointer<vtkDataArray> geoArray = dataSet->GetPointData()->GetArray("geo");

   int numberOfPoints = dataSet->GetNumberOfPoints();

   double x[3];
   int   ix[3];

   print("create matrix: start");
   timer_create->StartTimer();

   for (int i = 0; i < numberOfPoints; i++)
   {
      dataSet->GetPoint(i, x);
      get_node_indexes(x, ix);
      VMATRIX[index(Vx, ix[0], ix[1], ix[2])] = mvxArray->GetTuple1(i);
      VMATRIX[index(Vxx, ix[0], ix[1], ix[2])] = vxxArray->GetTuple1(i);
      VMATRIX[index(Vyy, ix[0], ix[1], ix[2])] = vyyArray->GetTuple1(i);
      VMATRIX[index(Vzz, ix[0], ix[1], ix[2])] = vzzArray->GetTuple1(i);
      FMATRIX[index2(ix[0], ix[1], ix[2])] = (int)fArray->GetTuple1(i);
      GMATRIX[index2(ix[0], ix[1], ix[2])] = (int)geoArray->GetTuple1(i);
   }

   print("create matrix: end");
   timer_create->StopTimer();
   print("create matrix time: "+toString(timer_create->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void read_ZDATA(string input, int nx)
{
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   print("read z data to "+input+": start");
   timer_read->StartTimer();
   std::ifstream istr;

   NX4 = nx;

   vector<double>zdatac(nx, 0.0);
   vector<double>zdatavx(nx, 0.0);
   vector<double>zdatavxx(nx, 0.0);
   vector<double>zdatavyy(nx, 0.0);
   vector<double>zdatavzz(nx, 0.0);

   istr.open(input.c_str(), std::ifstream::in);
   if(!istr)
   { 
      istr.clear();
   }
   string dummy;

   istr >> dummy >> dummy; 
   istr >> dummy >> dummy; 
   istr >> dummy; 

   int i = 0;

   while(istr.good())
   {
      istr >> zdatac[i];
      istr >> zdatavx[i] ;
      istr >> zdatavxx[i];
      istr >> zdatavyy[i];
      istr >> zdatavzz[i];
      i++;
   }

   istr.close();

   vector<int> bad_index(NX4, 0);

   int numberOfPoints = (int)zdatac.size();

   for(int x3 = 0; x3 < numberOfPoints; x3++)
   {
      if (zdatavx[x3] != 0)
      {
         ZDATAC.push_back(zdatac[x3]); 
         ZDATAVX.push_back(zdatavx[x3]); 
         ZDATAVXX.push_back(zdatavxx[x3]);  
         ZDATAVYY.push_back(zdatavyy[x3]);  
         ZDATAVZZ.push_back(zdatavzz[x3]);  
      }
   }

   NX4 = (int)ZDATAVX.size();

   print("read data set: end");
   timer_read->StopTimer();
   print("read z data time: "+toString(timer_read->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void v_averaging_z(double L)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   print("volume averaging: start");
   timer_averaging->StartTimer();

   print("L = "+toString(L));

   int l = cint(L/DELTAX);

   print("l = "+toString(l));

   double mdv;

   int numberOfPoints = (int)ZDATAVX.size();

   for(int x3 = 0; x3 < numberOfPoints; x3++)
   {
      double mvx = 0;
      double mvxx = 0;
      double mvyy = 0;
      double mvzz = 0;

      for(int z=-l;z<=+l;z++)
      {
         int xx = 0;
         int yy = 0;
         int zz = x3+z;
         correct_index(xx, yy,zz);

         mdv = m2((double)z, (double)l);
         mvx  += mdv*ZDATAVX[zz] ;
         mvxx += mdv*ZDATAVXX[zz];
         mvyy += mdv*ZDATAVYY[zz];
         mvzz += mdv*ZDATAVZZ[zz];

      }

      ZDATAVX[x3]  = mvx ;
      ZDATAVXX[x3] = mvxx;
      ZDATAVYY[x3] = mvyy;
      ZDATAVZZ[x3] = mvzz;
   }

timer_averaging->StopTimer();
print("volume averaging: end");
print("volume averaging time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void VolumeAveragingWithMatrix(std::string dataName, std::string dataNameL, std::string dataNameRMS, int numberOfLevels, /*int msize[3],*/ double L, double deltax, std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
   timer_global->StartTimer();
   //DELTAX = deltax;

   DELTAX = deltax / (double)(1<<(numberOfLevels-1));

   //create_matrix(dataName, dataNameL, dataNameRMS, numberOfLevels);
   //write_matrix_to_image_file(output+"s");
   //xy_averaging();
   //v_averaging2(L);
   //write_matrix_to_image_file(output+"v");

   //string file = output + ".vti";
   //read_image_data(file);

   //s_averaging();
   //write_ZDATA(output);
   //read_ZDATA(output+".dat", 68);

   read_ZDATA("d:/Projects/SFB880/BKanal/newData/Ludwig/BKanalAllVolXY.dat", 604);

   v_averaging_z(L);

   write_ZDATA(output+"ZDATA");

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}



