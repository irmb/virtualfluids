#ifndef FEADHOCTRIFACEMESH3D_H
#define FEADHOCTRIFACEMESH3D_H

#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

//extern "C"
//{
   //#include "mshpi.h"
    #include "fsi_interface.h"
    #include "fsi_user_interface.h"
//}

#include "./FeTriFaceMesh3D.h"

#ifdef CAB_RCF
#include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF


/*=========================================================================*/
/* FeAdhocTriFaceMesh3D                                                    */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
*/
class FeAdhocTriFaceMesh3D : public FeTriFaceMesh3D
{
public:
   FeAdhocTriFaceMesh3D();
   FeAdhocTriFaceMesh3D(std::string name, Mesh *mesh);

   Mesh*  getMesh() { return mesh; }

   void writeAdhoCMeshForStefan(string filename);
   std::vector< ::Vertex*>* getAdhocVertices() { return this->adhocVertices; }

#ifdef CAB_RCF
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      serializeParent<FeTriFaceMesh3D>(ar, *this);
   }
#endif //CAB_RCF

private:
   Mesh* mesh;
   std::vector< ::Vertex*>* adhocVertices;
};

#ifdef RCF_USE_SF_SERIALIZATION
UB_AUTO_RUN_NAMED(   SF::registerType<FeAdhocTriFaceMesh3D  >("FeAdhocTriFaceMesh3D  ")    , SF_FeAdhocTriFaceMesh3D     );
UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< FeTriFaceMesh3D, FeAdhocTriFaceMesh3D>() ), SF_FeAdhocTriFaceMesh3D_BD1 );
UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbTriFaceMesh3D, FeAdhocTriFaceMesh3D>() ), SF_FeAdhocTriFaceMesh3D_BD2 );
UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, FeAdhocTriFaceMesh3D>() ), SF_FeAdhocTriFaceMesh3D_BD3 );
#endif //RCF_USE_SF_SERIALIZATION

#endif //FEADHOCTRIFACEMESH3D_H
