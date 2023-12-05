#ifndef WRITEVTKSG_HPP
#define WRITEVTKSG_HPP

#include <stdio.h>
#include <fstream>
#include <sstream>
// #include <math.h>
#include <cmath>
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"
#include "basics/utilities/UbSystem.h"

//using namespace std;

namespace VtkSGWriter
{

void writeVTKsg(unsigned int nx, unsigned int ny, unsigned int nz, int startoffx, int startoffy, int startoffz,
                unsigned int nnx, unsigned int nny, unsigned int nnz, unsigned int start_z, unsigned int end_z, std::string& fname,
                   unsigned int*   bcMatH, float* vxH, float* vyH, float* vzH, float* rhoH, float v_ratio, float rho_ratio, 
                float startx, float starty, float startz, float deltax)
{


    unsigned int i,j,k,m;//,x,y,z;
   //unsigned int l;
    unsigned int nn;
    unsigned int tmp;
   float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
   //unsigned int tmp1;
   //float rho,ux,uy,uz,m6;
    //float f[19];

    float uzMax=0.0;

    nn = nnx*nny*nnz;

    ostringstream ostr;
    ofstream ofile(fname.c_str());

   ostr << "# vtk DataFile Version 3.0\n";
   ostr << fname << "\n";
   ostr << "ASCII\n";
   ostr << "DATASET STRUCTURED_POINTS\n";
   ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << "\n";
   ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
   ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

   ostr << "POINT_DATA " << nn << "\n";

   ostr << "SCALARS Geo float\n";
   ostr << "LOOKUP_TABLE default\n";

   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            tmp = bcMatH[m];
            ostr << tmp << "\n";
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   ostr << "SCALARS Press float\n";
   ostr << "LOOKUP_TABLE default\n";
   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            if( bcMatH[m] == GEO_SOLID){
               tmp_rho=0.0f;//1.0f/3.0f;
            }
            else{
               tmp_rho=rhoH[m] / 3.0f * rho_ratio * v_ratio *v_ratio;

            }
            ostr << tmp_rho << "\n";
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   ostr << "VECTORS Speed float\n";
   //ostr << "LOOKUP_TABLE default\n";
   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            if( bcMatH[m] == GEO_SOLID){
               tmp_ux=0.0;
               tmp_uy=0.0;
               tmp_uz=0.0;
            }
            else{
               tmp_ux=vxH[m] * v_ratio;
               tmp_uy=vyH[m] * v_ratio;
               tmp_uz=vzH[m] * v_ratio;
            }
            ostr << tmp_ux << " " << tmp_uy << " " << tmp_uz << "\n";
            if( std::fabs(tmp_uz) > std::fabs(uzMax) ) uzMax = tmp_uz;
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

    printf("uzMax: %.8f", uzMax);
}



void writeVTKsgSP(unsigned int nx_D, unsigned int ny_D, unsigned int nz_D, int startoffx_D, int startoffy_D, int startoffz_D,
                  unsigned int nnx_D, unsigned int nny_D, unsigned int nnz_D, unsigned int start_z, unsigned int end_z, std::string& fname,
                  unsigned int*   bcMatH, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real* pressH_D, real v_ratio_D, real rho_ratio_D, 
                  real startx_D, real starty_D, real startz_D, real deltax_D, unsigned int *kFull)
 {
    //Test
    real* /*float* double**/ vxH      = vxH_D      ;
    real* /*float* double**/ vyH      = vyH_D      ;
    real* /*float* double**/ vzH      = vzH_D      ;
    real* /*float* double**/ rhoH     = rhoH_D     ;
    real* /*float* double**/ pressH   = pressH_D   ;
    double  v_ratio  = v_ratio_D  ;
    double  rho_ratio= rho_ratio_D;
    double  startx   = startx_D   ;
    double  starty   = starty_D   ; 
    double  startz   = startz_D   ;//(double)start_z;
    double  deltax   = deltax_D   ;
    int startoffx = startoffx_D;
    int startoffy = startoffy_D;
    int startoffz = startoffz_D;
    unsigned int nnx = nnx_D;
    unsigned int nny = nny_D;
    unsigned int nnz = nnz_D;
    unsigned int nx = nx_D;
    unsigned int ny = ny_D;
    unsigned int nz = nz_D;
    ////f�r Geller
    //startoffx = 0;
    //startoffy = 0;
    //startoffz = 0;
    //nx = nnx;
    //ny = nny;
    //nz = nnz;
    //deltax = 410.0/48.0;

    unsigned int i,j,k,m;
    unsigned int nn;
    unsigned int tmp;
    double tmp_rho,tmp_ux,tmp_uy,tmp_uz;


    double uzMax=0.0;

    nn = nnx*nny*nnz;

    ostringstream ostr;
    ofstream ofile(fname.c_str());

    ostr << "# vtk DataFile Version 3.0\n";
    ostr << fname << "\n";
    ostr << "ASCII\n";
    ostr << "DATASET STRUCTURED_POINTS\n";
    ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << "\n";
    ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
    ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

    ostr << "POINT_DATA " << nn << "\n";
    //ostr << "CELL_DATA " << nn << "\n";

    ostr << "SCALARS Geo int\n";
    ostr << "LOOKUP_TABLE default\n";

    for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
        for(j=startoffy ; j<=nny+startoffy-1 ; j++){
            for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                m = nx*(ny*k + j) + i;
                if(kFull[m]==0) tmp = 99;
                //else tmp = bcMatH[kFull[m]-1];//achtung!!! -1 
                tmp = bcMatH[kFull[m]];
                //tmp = bcMatH[m];
                ostr << tmp << "\n";
            }
        }
    }
    ofile<<ostr.str();
    ostr.str("");
    ostr.clear();

    ostr << "SCALARS Press double\n";
    ostr << "LOOKUP_TABLE default\n";
    for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
       for(j=startoffy ; j<=nny+startoffy-1 ; j++){
          for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
             m = nx*(ny*k + j) + i;
             if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                tmp_rho=0.0;//1.0f/3.0f;
             }
             else{
                tmp_rho=rhoH[kFull[m]] / 3.0f * rho_ratio * v_ratio *v_ratio;

             }
             ostr << tmp_rho << "\n";
          }
       }
    }
    ofile<<ostr.str();
    ostr.str("");
    ostr.clear();

    ostr << "SCALARS PressM double\n";
    ostr << "LOOKUP_TABLE default\n";
    for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
       for(j=startoffy ; j<=nny+startoffy-1 ; j++){
          for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
             m = nx*(ny*k + j) + i;
             if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                tmp_rho=0.0;//1.0f/3.0f;
             }
             else{
                tmp_rho=pressH[kFull[m]] / 3.0f * rho_ratio * v_ratio *v_ratio;

             }
             ostr << tmp_rho << "\n";
          }
       }
    }
    ofile<<ostr.str();
    ostr.str("");
    ostr.clear();

    ostr << "VECTORS Speed double\n";
    for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
       for(j=startoffy ; j<=nny+startoffy-1 ; j++){
          for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
             m = nx*(ny*k + j) + i;
             if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                tmp_ux=0.0;
                tmp_uy=0.0;
                tmp_uz=0.0;
             }
             else{
                tmp_ux=vxH[kFull[m]] * v_ratio;
                tmp_uy=vyH[kFull[m]] * v_ratio;
                tmp_uz=vzH[kFull[m]] * v_ratio;
             }
             ostr << tmp_ux << " " << tmp_uy << " " << tmp_uz << "\n";
             if( std::fabs(tmp_uz) > std::fabs(uzMax) ) uzMax = tmp_uz;
          }
       }
    }
    ofile<<ostr.str();
    ostr.str("");
    ostr.clear();

    printf(" \n uzMax: %.8f", (float)uzMax);
 }








 void writeVTKsgSPbin(unsigned int nx_D, unsigned int ny_D, unsigned int nz_D, int startoffx_D, int startoffy_D, int startoffz_D,
                     unsigned int nnx_D, unsigned int nny_D, unsigned int nnz_D, unsigned int start_z, unsigned int end_z, std::string& fname,
                     unsigned int*   bcMatH, int*   bcMatHFull, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real* pressH_D, real v_ratio_D, real rho_ratio_D, 
                     real startx_D, real starty_D, real startz_D, real deltax_D, unsigned int *kFull)
 {
     //Test 
     real* /*float* double**/ vxH      = vxH_D      ;
     real* /*float* double**/ vyH      = vyH_D      ;
     real* /*float* double**/ vzH      = vzH_D      ;
     real* /*float* double**/ rhoH     = rhoH_D     ;
     real* /*float* double**/ pressH   = pressH_D   ;
     float  v_ratio  = (float)v_ratio_D  ;
     float  rho_ratio= (float)rho_ratio_D;
     float  startx   = (float)startx_D   ;
     float  starty   = (float)starty_D   ; 
     float  startz   = (float)startz_D   ;//(double)start_z;
     float  deltax   = (float)deltax_D   ;
     int startoffx = startoffx_D;
     int startoffy = startoffy_D;
     int startoffz = startoffz_D;
     unsigned int nnx = nnx_D;
     unsigned int nny = nny_D;
     unsigned int nnz = nnz_D;
     unsigned int nx = nx_D;
     unsigned int ny = ny_D;
     unsigned int nz = nz_D;
     unsigned int i,j,k,m;
     unsigned int nn;
     //unsigned int tmp;
     //double tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     int tmp;
     bool swapByte = UbSystem::isLittleEndian();
     float uxMax=0.0f;


     nn = nnx*nny*nnz;

     ofstream ostr(fname.c_str(),ios::out | ios::binary);

     ostr << "# vtk DataFile Version 3.0\n";
     ostr << fname << "\n";
     ostr << "BINARY\n";
     ostr << "DATASET STRUCTURED_POINTS\n";
     ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << " \n";
     ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
     ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

     ostr << "POINT_DATA " << nn << "\n";

     ostr << "SCALARS Geo int\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 //if(kFull[m]==0) tmp = 16;
                 tmp = bcMatH[kFull[m]];
                 //tmp = bcMatHFull[m];
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp,sizeof(int));
                 ostr.write((char*)&tmp, sizeof(int));
             }
         }
     }

     ostr << "\n";
     ostr << "SCALARS Press float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_rho=0.0;//1.0f/3.0f;
                 }
                 else{
                     tmp_rho=(float)(rhoH[kFull[m]] / 3.0f * rho_ratio * v_ratio *v_ratio);
                     //tmp_rho=rhoH[kFull[m]] / 3.0f ;

                 }
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\n";
     ostr << "SCALARS PressM float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_rho=0.0f;//1.0f/3.0f;
                 }
                 else{
                     tmp_rho=(float)(pressH[kFull[m]] / 3.0f * rho_ratio * v_ratio *v_ratio);
                     //tmp_rho=pressH[kFull[m]] / 3.0f ;

                 }
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\n";
     ostr << "VECTORS Speed float\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_ux=0.0;
                     tmp_uy=0.0;
                     tmp_uz=0.0;
                 }
                 else{
                     tmp_ux=(float)(vxH[kFull[m]] * v_ratio);
                     tmp_uy=(float)(vyH[kFull[m]] * v_ratio);
                     tmp_uz=(float)(vzH[kFull[m]] * v_ratio);
                 }
                 if( std::fabs(tmp_ux) > std::fabs(uxMax) ) uxMax = tmp_ux;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_ux,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uy,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uz,sizeof(float));
                 }
                 ostr.write((char*)&tmp_ux, sizeof(float));
                 ostr.write((char*)&tmp_uy, sizeof(float));
                 ostr.write((char*)&tmp_uz, sizeof(float));
             }
         }
     }
     printf(" \n uxMax: %.8f", (float)uxMax);
 }



void writeVTKmedSPbin(   unsigned int nx_D, unsigned int ny_D, unsigned int nz_D, int startoffx_D, int startoffy_D, int startoffz_D,
                         unsigned int nnx_D, unsigned int nny_D, unsigned int nnz_D, unsigned int start_z, unsigned int end_z, std::string& fname,
                         unsigned int*   bcMatH, int*   bcMatHFull, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real* pressH_D, real tdiff, real v_ratio_D, real rho_ratio_D, 
                         real startx_D, real starty_D, real startz_D, real deltax_D, unsigned int *kFull)
 {
     //Test 
     real* /*float* double**/ vxH      = vxH_D      ;
     real* /*float* double**/ vyH      = vyH_D      ;
     real* /*float* double**/ vzH      = vzH_D      ;
     real* /*float* double**/ rhoH     = rhoH_D     ;
     real* /*float* double**/ pressH   = pressH_D   ;
     float  v_ratio  = (float)v_ratio_D  ;
     float  rho_ratio= (float)rho_ratio_D;
     float  startx   = (float)startx_D   ;
     float  starty   = (float)starty_D   ; 
     float  startz   = (float)startz_D   ;//(double)start_z;
     float  deltax   = (float)deltax_D   ;
     int startoffx = startoffx_D;
     int startoffy = startoffy_D;
     int startoffz = startoffz_D;
     unsigned int nnx = nnx_D;
     unsigned int nny = nny_D;
     unsigned int nnz = nnz_D;
     unsigned int nx = nx_D;
     unsigned int ny = ny_D;
     unsigned int nz = nz_D;
     unsigned int i,j,k,m;
     unsigned int nn;
     //unsigned int tmp;
     //double tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     int tmp;
     bool swapByte = UbSystem::isLittleEndian();
     float uxMax=0.0f;


     nn = nnx*nny*nnz;

     ofstream ostr(fname.c_str(),ios::out | ios::binary);

     ostr << "# vtk DataFile Version 3.0\n";
     ostr << fname << "\n";
     ostr << "BINARY\n";
     ostr << "DATASET STRUCTURED_POINTS\n";
     ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << " \n";
     ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
     ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

     ostr << "POINT_DATA " << nn << "\n";

     ostr << "SCALARS Geo int\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 //if(kFull[m]==0) tmp = 16;
                 tmp = bcMatH[kFull[m]];
                 //tmp = bcMatHFull[m];
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp,sizeof(int));
                 ostr.write((char*)&tmp, sizeof(int));
             }
         }
     }

     ostr << "\n";
     ostr << "SCALARS Press float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_rho=0.0;//1.0f/3.0f;
                 }
                 else{
                     tmp_rho=(float)((rhoH[kFull[m]] / tdiff) / 3.0f * rho_ratio * v_ratio *v_ratio );
                     //tmp_rho=rhoH[kFull[m]] / 3.0f ;

                 }
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\n";
     ostr << "SCALARS PressM float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_rho=0.0f;//1.0f/3.0f;
                 }
                 else{
                     tmp_rho=(float)((pressH[kFull[m]] / tdiff) / 3.0f * rho_ratio * v_ratio *v_ratio );
                     //tmp_rho=pressH[kFull[m]] / 3.0f ;

                 }
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\n";
     ostr << "VECTORS Speed float\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 m = nx*(ny*k + j) + i;
                 if( (bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
                     tmp_ux=0.0;
                     tmp_uy=0.0;
                     tmp_uz=0.0;
                 }
                 else{
                     tmp_ux=(float)((vxH[kFull[m]] / tdiff) * v_ratio);
                     tmp_uy=(float)((vyH[kFull[m]] / tdiff) * v_ratio);
                     tmp_uz=(float)((vzH[kFull[m]] / tdiff) * v_ratio);
                 }
                 if( std::fabs(tmp_ux) > std::fabs(uxMax) ) uxMax = tmp_ux;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_ux,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uy,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uz,sizeof(float));
                 }
                 ostr.write((char*)&tmp_ux, sizeof(float));
                 ostr.write((char*)&tmp_uy, sizeof(float));
                 ostr.write((char*)&tmp_uz, sizeof(float));
             }
         }
     }
     printf(" \n uxMax: %.8f", (float)uxMax);
 }



 void writeVTKsgSPbinTEST(std::string& fname)
 {
     unsigned int i,j,k;
     float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     int tmp;
     bool swapByte = UbSystem::isLittleEndian();


     ofstream ostr(fname.c_str(),ios::out | ios::binary);

     ostr << "# vtk DataFile Version 3.0\n";
     ostr << fname << "\n";
     ostr << "BINARY\n";
     ostr << "DATASET STRUCTURED_POINTS\n";
     ostr << "DIMENSIONS " << 10 << " " << 10 << " " << 10 << " \n";
     ostr << "ORIGIN " << 0 << " " << 0 << " " << 0 << " \n";
     ostr << "SPACING " << 1 << " " << 1 << " " << 1 <<" \n";

     ostr << "POINT_DATA " << 1000 << "\n";

     ostr << "SCALARS Geo int\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=0 ; k<10 ; k++){
         for(j=0 ; j<10 ; j++){
             for(i=0 ; i<10 ; i++){
                 tmp = 1;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp,sizeof(int));
                 }
                 ostr.write((char*)&tmp, sizeof(int));
             }
         }
     }

     ostr << "\nSCALARS Press float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=0 ; k<10 ; k++){
         for(j=0 ; j<10 ; j++){
             for(i=0 ; i<10 ; i++){
                 tmp_rho=0.1f;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 }
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\nSCALARS PressM float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=0 ; k<10 ; k++){
         for(j=0 ; j<10 ; j++){
             for(i=0 ; i<10 ; i++){
                 tmp_rho=0.2f;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 }
                 ostr.write((char*)&tmp_rho, sizeof(float));
             }
         }
     }

     ostr << "\nVECTORS Speed float\n";
     for(k=0 ; k<10 ; k++){
         for(j=0 ; j<10 ; j++){
             for(i=0 ; i<10 ; i++){
                 tmp_ux=0.01f;
                 tmp_uy=0.02f;
                 tmp_uz=0.03f;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_ux,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uy,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uz,sizeof(float));
                 }
                 ostr.write((char*)&tmp_ux, sizeof(float));
                 ostr.write((char*)&tmp_uy, sizeof(float));
                 ostr.write((char*)&tmp_uz, sizeof(float));
             }
         }
     }
 }







void writeVTKsgThS(unsigned int nx, unsigned int ny, unsigned int nz, int startoffx, int startoffy, int startoffz,
                   unsigned int nnx, unsigned int nny, unsigned int nnz, unsigned int start_z, unsigned int end_z, std::string& fname,
                   unsigned int*   bcMatH, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real v_ratio_D, real rho_ratio_D, 
                   real startx_D, real starty_D, real startz_D, real deltax_D, unsigned int *kFull, real* ConcD)
{
   //Test
    real* /*float* double**/ vxH      = vxH_D      ;
    real* /*float* double**/ vyH      = vyH_D      ;
    real* /*float* double**/ vzH      = vzH_D      ;
    real* /*float* double**/ rhoH     = rhoH_D     ;
    real* /*float* double**/ ConcH    = ConcD   ;
   double  v_ratio  = v_ratio_D  ;
   double  rho_ratio= rho_ratio_D;
   double  startx   = startx_D   ;
   double  starty   = starty_D   ; 
   double  startz   = startz_D   ;//(double)start_z;
   double  deltax   = deltax_D   ;

   unsigned int i,j,k,m;
   unsigned int nn;
   unsigned int tmp;
   double tmp_rho,tmp_ux,tmp_uy,tmp_uz, temp_Conc;


   double uzMax=0.0;

   nn = nnx*nny*nnz;

   ostringstream ostr;
   ofstream ofile(fname.c_str());

   ostr << "# vtk DataFile Version 3.0\n";
   ostr << fname << "\n";
   ostr << "ASCII\n";
   ostr << "DATASET STRUCTURED_POINTS\n";
   ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << "\n";
   ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
   ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

   ostr << "POINT_DATA " << nn << "\n";

   ostr << "SCALARS Geo int\n";
   ostr << "LOOKUP_TABLE default\n";

   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            tmp = bcMatH[kFull[m]];
            ostr << tmp << "\n";
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   ostr << "SCALARS Press float\n";
   ostr << "LOOKUP_TABLE default\n";
   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
               tmp_rho=0.0;
            }
            else{
               tmp_rho=rhoH[kFull[m]] / 3.0f * rho_ratio * v_ratio *v_ratio;

            }
            ostr << tmp_rho << "\n";
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   ostr << "SCALARS Conc float\n";
   ostr << "LOOKUP_TABLE default\n";
   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
               temp_Conc=0.0;
            }
            else{
               temp_Conc=ConcH[kFull[m]];

            }
            ostr << temp_Conc << "\n";
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   ostr << "VECTORS Speed float\n";
   for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
      for(j=startoffy ; j<=nny+startoffy-1 ; j++){
         for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
            m = nx*(ny*k + j) + i;
            if((bcMatH[kFull[m]] == GEO_SOLID) && (bcMatH[kFull[m]] == GEO_VOID)){
               tmp_ux=0.0;
               tmp_uy=0.0;
               tmp_uz=0.0;
            }
            else{
               tmp_ux=vxH[kFull[m]] * v_ratio;
               tmp_uy=vyH[kFull[m]] * v_ratio;
               tmp_uz=vzH[kFull[m]] * v_ratio;
            }
            ostr << tmp_ux << " " << tmp_uy << " " << tmp_uz << "\n";
            if( std::fabs(tmp_uz) > std::fabs(uzMax) ) uzMax = tmp_uz;
         }
      }
   }
   ofile<<ostr.str();
   ostr.str("");
   ostr.clear();

   printf(" \n uzMax: %.8f", (float)uzMax);
 }





 void writeVTKsgSPbinAS(unsigned int nx_D, unsigned int ny_D, unsigned int nz_D, int startoffx_D, int startoffy_D, int startoffz_D,
                        unsigned int nnx_D, unsigned int nny_D, unsigned int nnz_D, unsigned int start_z, unsigned int end_z, std::string& fname,
                        unsigned int*   bcMatH, int*   bcMatHFull, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real* pressH_D, 
                        real v_ratio_D, real rho_ratio_D, 
                        real startx_D, real starty_D, real startz_D, real deltax_D, 
                        real *coordX, real *coordY, real *coordZ)
 {
     //Test 
     real* /*float* double**/ vxH      = vxH_D      ;
     real* /*float* double**/ vyH      = vyH_D      ;
     real* /*float* double**/ vzH      = vzH_D      ;
     real* /*float* double**/ rhoH     = rhoH_D     ;
     real* /*float* double**/ pressH   = pressH_D   ;
     float  v_ratio  = (float)v_ratio_D  ;
     float  rho_ratio= (float)rho_ratio_D;
     float  startx   = (float)startx_D   ;
     float  starty   = (float)starty_D   ; 
     float  startz   = (float)startz_D   ;//(double)start_z;
     float  deltax   = (float)deltax_D   ;
     int startoffx = startoffx_D;
     int startoffy = startoffy_D;
     int startoffz = startoffz_D;
     unsigned int nnx = nnx_D;
     unsigned int nny = nny_D;
     unsigned int nnz = nnz_D;
     unsigned int nx = nx_D;
     unsigned int ny = ny_D;
     unsigned int nz = nz_D;
     unsigned int i,j,k;//,m;
     unsigned int nn;
     float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     int tmp;
     bool swapByte = UbSystem::isLittleEndian();
     float uxMax=0.0f;

     unsigned int count=1;

     nn = nnx*nny*nnz;

     ofstream ostr(fname.c_str(),ios::out | ios::binary);

     ostr << "# vtk DataFile Version 3.0\n";
     ostr << fname << "\n";
     ostr << "BINARY\n";
     ostr << "DATASET STRUCTURED_POINTS\n";
     ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << " \n";
     ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
     ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

     ostr << "POINT_DATA " << nn << "\n";

     ostr << "SCALARS Geo int\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp = bcMatH[count]; 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp = 99;//GEO_VOID; 
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp,sizeof(int));
                 ostr.write((char*)&tmp, sizeof(int));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "SCALARS Press float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_rho = (float)(rhoH[count] / 3.0f * rho_ratio * v_ratio *v_ratio); 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp_rho=0.0f;
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "SCALARS PressM float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_rho = (float)(pressH[count] / 3.0f * rho_ratio * v_ratio *v_ratio); 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp_rho=0.0f;
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "VECTORS Speed float\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_ux=(float)(vxH[count] * v_ratio);
                     tmp_uy=(float)(vyH[count] * v_ratio);
                     tmp_uz=(float)(vzH[count] * v_ratio);
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else {
                     tmp_ux=0.0f;
                     tmp_uy=0.0f;
                     tmp_uz=0.0f;
                 }
                 //////////////////////////////////////////////////////////////////////////
                 if( std::fabs(tmp_ux) > std::fabs(uxMax) ) uxMax = tmp_ux;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_ux,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uy,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uz,sizeof(float));
                 }
                 ostr.write((char*)&tmp_ux, sizeof(float));
                 ostr.write((char*)&tmp_uy, sizeof(float));
                 ostr.write((char*)&tmp_uz, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }
     printf(" \n uxMax: %.8f", (float)uxMax);
 }



 void writeVTKmedSPbinAS(unsigned int nx_D, unsigned int ny_D, unsigned int nz_D, int startoffx_D, int startoffy_D, int startoffz_D,
                         unsigned int nnx_D, unsigned int nny_D, unsigned int nnz_D, unsigned int start_z, unsigned int end_z, std::string& fname,
                         unsigned int*   bcMatH, int*   bcMatHFull, real* vxH_D, real* vyH_D, real* vzH_D, real* rhoH_D, real* pressH_D, 
                         real tdiff, real v_ratio_D, real rho_ratio_D, 
                         real startx_D, real starty_D, real startz_D, real deltax_D, 
                         unsigned int *coordX, unsigned int *coordY, unsigned int *coordZ)
 {
     //Test 
     real* /*float* double**/ vxH      = vxH_D      ;
     real* /*float* double**/ vyH      = vyH_D      ;
     real* /*float* double**/ vzH      = vzH_D      ;
     real* /*float* double**/ rhoH     = rhoH_D     ;
     real* /*float* double**/ pressH   = pressH_D   ;
     float  v_ratio  = (float)v_ratio_D  ;
     float  rho_ratio= (float)rho_ratio_D;
     float  startx   = (float)startx_D   ;
     float  starty   = (float)starty_D   ; 
     float  startz   = (float)startz_D   ;//(double)start_z;
     float  deltax   = (float)deltax_D   ;
     int startoffx = startoffx_D;
     int startoffy = startoffy_D;
     int startoffz = startoffz_D;
     unsigned int nnx = nnx_D;
     unsigned int nny = nny_D;
     unsigned int nnz = nnz_D;
     unsigned int nx = nx_D;
     unsigned int ny = ny_D;
     unsigned int nz = nz_D;
     unsigned int i,j,k;//,m;
     unsigned int nn;
     //unsigned int tmp;
     //double tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     float tmp_rho,tmp_ux,tmp_uy,tmp_uz;
     int tmp;
     bool swapByte = UbSystem::isLittleEndian();
     float uxMax=0.0f;

     unsigned int count=1;

     nn = nnx*nny*nnz;

     ofstream ostr(fname.c_str(),ios::out | ios::binary);

     ostr << "# vtk DataFile Version 3.0\n";
     ostr << fname << "\n";
     ostr << "BINARY\n";
     ostr << "DATASET STRUCTURED_POINTS\n";
     ostr << "DIMENSIONS " << nnx << " " << nny << " " << nnz << " \n";
     ostr << "ORIGIN " << startx << " " << starty << " " << startz << " \n";
     ostr << "SPACING " << deltax << " " << deltax << " " << deltax <<" \n";

     ostr << "POINT_DATA " << nn << "\n";

     ostr << "SCALARS Geo int\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp = bcMatH[count]; 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp = GEO_VOID; 
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp,sizeof(int));
                 ostr.write((char*)&tmp, sizeof(int));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "SCALARS Press float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_rho = (float)((rhoH[count] / tdiff) / 3.0f * rho_ratio * v_ratio *v_ratio); 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp_rho=0.0f;
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "SCALARS PressM float\n";
     ostr << "LOOKUP_TABLE default\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_rho = (float)((pressH[count] / tdiff) / 3.0f * rho_ratio * v_ratio *v_ratio); 
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else tmp_rho=0.0f;
                 //////////////////////////////////////////////////////////////////////////
                 if(swapByte) UbSystem::swapByteOrder((unsigned char*)&tmp_rho,sizeof(float));
                 ostr.write((char*)&tmp_rho, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }

     //reset var "count"
     count=1;

     ostr << "\n";
     ostr << "VECTORS Speed float\n";
     for(k=startoffz ; k<=nnz+startoffz-1 ; k++){
         for(j=startoffy ; j<=nny+startoffy-1 ; j++){
             for(i=startoffx ; i<=nnx+startoffx-1 ; i++){
                 //////////////////////////////////////////////////////////////////////////
                 if(coordX[count]==i && coordY[count]==j && coordZ[count]==k) {
                     tmp_ux=(float)((vxH[count] / tdiff) * v_ratio);
                     tmp_uy=(float)((vyH[count] / tdiff) * v_ratio);
                     tmp_uz=(float)((vzH[count] / tdiff) * v_ratio);
                     count++; //count wird nur im "if" erh�ht. stellt sicher, dass wir im array nur einen schritt weiter gehen wenn eine kombination gefunden wurde
                 }
                 else {
                     tmp_ux=0.0f;
                     tmp_uy=0.0f;
                     tmp_uz=0.0f;
                 }
                 //////////////////////////////////////////////////////////////////////////
                 if( std::fabs(tmp_ux) > std::fabs(uxMax) ) uxMax = tmp_ux;
                 if(swapByte)
                 {
                     UbSystem::swapByteOrder((unsigned char*)&tmp_ux,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uy,sizeof(float));
                     UbSystem::swapByteOrder((unsigned char*)&tmp_uz,sizeof(float));
                 }
                 ostr.write((char*)&tmp_ux, sizeof(float));
                 ostr.write((char*)&tmp_uy, sizeof(float));
                 ostr.write((char*)&tmp_uz, sizeof(float));
                 //////////////////////////////////////////////////////////////////////////
             }
         }
     }
     printf(" \n uxMax: %.8f", (float)uxMax);
 }




}

#endif
