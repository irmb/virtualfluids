#ifndef FETRIFACEMESH3D_H
#define FETRIFACEMESH3D_H

#include <sstream>
#include <iostream>
#include <vector>
#include <map>

#include "./../GbTriFaceMesh3D.h"

#ifdef CAB_RCF
#include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

/*=========================================================================*/
/* FeTriFaceMesh3D                                                                  */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
*/
class FeTriFaceMesh3D : public GbTriFaceMesh3D
{                             
public:
   class VertexAttributes
   {
   private:
      double Fx;
      double Fy;
      double Fz;
      double sumFx;
      double sumFy;
      double sumFz;
      double velocityX;
      double velocityY;
      double velocityZ;
      double normalX;
      double normalY;
      double normalZ;
      UbTupleDouble6 stresses; 
      double area;
      double p; //pressure
   public:
      VertexAttributes()
      {
         this->init();
      }
      ~VertexAttributes() {}

#ifdef CAB_RCF
      template<class Archive>
      void SF_SERIALIZE(Archive & ar)
      {
         ar & Fx; ar & Fy; ar & Fz;
         ar & velocityX; ar & velocityY; ar & velocityZ;
      }
#endif //CAB_RCF

      void init()
      {
         this->Fx         = 0.0;
         this->Fy         = 0.0;
         this->Fz         = 0.0;
         this->sumFx      = 0.0;
         this->sumFy      = 0.0;
         this->sumFz      = 0.0;
         val<1>(stresses) = 0.0;
         val<2>(stresses) = 0.0;
         val<3>(stresses) = 0.0;
         val<4>(stresses) = 0.0;
         val<5>(stresses) = 0.0;
         val<6>(stresses) = 0.0;
         this->area       = 0.0;
         this->p = 0.0;
         this->velocityX  = 0.0;
         this->velocityY  = 0.0;
         this->velocityZ  = 0.0;
         this->normalX  = 0.0;
         this->normalY  = 0.0;
         this->normalZ  = 0.0;
      }
      double getVelocityX()         { return this->velocityX; }
      double getVelocityY()         { return this->velocityY; }
      double getVelocityZ()         { return this->velocityZ; }
      void   setVelocityX(double x) { this->velocityX = x;    }
      void   setVelocityY(double y) { this->velocityY = y;    }
      void   setVelocityZ(double z) { this->velocityZ = z;    }

      double getNormalX()         { return this->normalX; }
      double getNormalY()         { return this->normalY; }
      double getNormalZ()         { return this->normalZ; }
      void   setNormalX(double x) { this->normalX = x;    }
      void   setNormalY(double y) { this->normalY = y;    }
      void   setNormalZ(double z) { this->normalZ = z;    }

      double getFX()          { return this->Fx; }
      double getFY()          { return this->Fy; }
      double getFZ()          { return this->Fz; }
      void   setFX(double FX) { this->Fx = FX;   }
      void   setFY(double FY) { this->Fy = FY;   }
      void   setFZ(double FZ) { this->Fz = FZ;   }
      void   addFX(double FX) { this->Fx += FX;  }
      void   addFY(double FY) { this->Fy += FY;  }
      void   addFZ(double FZ) { this->Fz += FZ;  }

      double getSumFX()          { return this->sumFx; }
      double getSumFY()          { return this->sumFy; }
      double getSumFZ()          { return this->sumFz; }
      void   setSumFX(double FX) { this->sumFx = FX;   }
      void   setSumFY(double FY) { this->sumFy = FY;   }
      void   setSumFZ(double FZ) { this->sumFz = FZ;   }
      void   addSumFX(double FX) { this->sumFx += FX;  }
      void   addSumFY(double FY) { this->sumFy += FY;  }
      void   addSumFZ(double FZ) { this->sumFz += FZ;  }

      UbTupleDouble6& getStresses() { return this->stresses; }
      
      double getArea()            { return this->area;  }
      void   setArea(double area) { this->area  = area; }
      void   addArea(double area) { this->area += area; }

      double getPressure()         { return this->p;  }
      void   setPressure(double p) { this->p = p; }

   };
/*=========================================================================*/
/*=========================================================================*/
/*=========================================================================*/
   //class VertexTriFaceMap : public std::multimap<Vertex*, TriFace*>
   //{
   //public:
   //   VertexTriFaceMap()  {}
   //   /*=========================================================================*/
   //   void setVertexTriFaceRelation(Vertex* v, TriFace* tri)
   //   {
   //      this->insert(std::pair<Vertex*,TriFace*>(v,tri));
   //   }
   //   /*=========================================================================*/
   //   int getNumberOfTriFaces(Vertex* v)
   //   {  
   //      return((int)this->count(v));
   //   }
   //   /*=========================================================================*/
   //   std::vector<TriFace*> getTriFacesForVertex(Vertex* v)
   //   {
   //      std::vector<TriFace*> trivector;
   //      unsigned number = (unsigned)this->count(v);
   //      std::multimap<Vertex*,TriFace*>::iterator mapIterator = this->find(v);
   //      for(unsigned u =0; u<number; u++) 
   //      {
   //         trivector.push_back(mapIterator->second);
   //         mapIterator ++;
   //      }
   //      return trivector;
   //   }
   //   //void deleteNeighbors(QtInteractor* interactor);
   //   ///*=========================================================================*/
   //};
/*=========================================================================*/
/*=========================================================================*/
/*=========================================================================*/
public:
   //#ifndef SWIG
   //VertexTriFaceMap vertexTriFaceMap;
   //#endif

   FeTriFaceMesh3D();
   FeTriFaceMesh3D(std::string name, std::vector<Vertex>* nodes, std::vector<TriFace>* triangles);

   std::vector<VertexAttributes>* getAttributes() { return this->attributes; }
   void resizeAttributes();
   //void createVertexTriFaceMap();

   UbTuple<std::string,std::string> writeMesh(std::string filename, WbWriter* writer, bool writeNormals=false, std::vector< std::string >* datanames=NULL, std::vector< std::vector < double > >* nodedata=NULL);

   static FeTriFaceMesh3D* createMeshByTriangles(std::string name, std::vector<GbTriangle3D*>* triangles);

   virtual ObObjectCreator* getCreator();

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbTriFaceMesh3D>(ar, *this);
      ar & attributes;
      //if( ArchiveTools::isReading(ar) ) this->createVertexTriFaceMap();
   }
#endif //CAB_RCF


protected:
   std::vector<VertexAttributes>* attributes;
   
};

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   UB_AUTO_RUN_NAMED(   SF::registerType<FeTriFaceMesh3D  >("FeTriFaceMesh3D  ")     , SF_FeTriFaceMesh3D     );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbTriFaceMesh3D, FeTriFaceMesh3D >() ), SF_FeTriFaceMesh3D_BD1 );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, FeTriFaceMesh3D >() ), SF_FeTriFaceMesh3D_BD2 );
#endif //RCF_USE_SF_SERIALIZATION


#endif //FETRIFACEMESH3D_H
