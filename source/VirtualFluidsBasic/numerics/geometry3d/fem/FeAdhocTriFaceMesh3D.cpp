#include "./FeAdhocTriFaceMesh3D.h"

#include "./../GbPoint3D.h"
#include "./../GbTriangle3D.h"
#include "./../GbTriangularMesh3D.h"

#include "./../../../basics/utilities/UbFileOutputASCII.h"
#include "./../../../basics/utilities/UbLogger.h"

/*====================================================*/
FeAdhocTriFaceMesh3D::FeAdhocTriFaceMesh3D():FeTriFaceMesh3D()
{
}
/*====================================================*/
FeAdhocTriFaceMesh3D::FeAdhocTriFaceMesh3D(std::string name, Mesh* mesh):FeTriFaceMesh3D()
{
   this->mesh = mesh;
   this->setName(name);

   std::cout << "FeAdhocTriFaceMesh3D:Konstruktor !!!"<<std::endl;
   std::cout << "num vertices: " << mesh->VL->nitem << ", num triangles: " << mesh->TL->nitem
             << ", num quads: " << mesh->QL->nitem << std::endl;

 //  this->writeAdhoCMeshForStefan("/scratch/geller/StudienMitAdhoC3D/mesh.inp");
   this->adhocVertices = new vector< ::Vertex*>;

   List *vertexlist  = mesh->VL;
   List *elementlist = mesh->TL;

   ListItem*  LI;
   ListItem*  LI2;
   Tri*       triangle;
 //  double z1, z2, z3;
   int    id1, id2, id3;
   ::Vertex *v1, *v2, *v3;

   //if (mesh->VL->status==open)  close_Vertex_List(mesh->VL);

   FOR_ALL(vertexlist->first, LI, 0)
   {
      ::Vertex* V = get_Vertex(LI);
      this->nodes->push_back(GbTriFaceMesh3D::Vertex((float)V->x, (float)V->y, (float)V->z));
      this->adhocVertices->push_back(V);
   }
   int countTris=0;
    int fred=0;
   FOR_ALL(elementlist->first, LI, 0)
   {
      triangle = get_Tri(LI);
    if(triangle==NULL) UBLOG(logINFO, "hugo - dreieck ist NULL");
      v1 = triangle->V[0];
      v2 = triangle->V[1];
      v3 = triangle->V[2];
      int count=0;
      id1=-1; id2=-1; id3=-1;
      FOR_ALL(vertexlist->first, LI2, 0)
      {
         ::Vertex* V = get_Vertex(LI2);
         if(v1==V) id1=count;
         if(v2==V) id2=count;
         if(v3==V) id3=count;
         if((id1!=-1) && (id2!=-1) && (id3!=-1))
         {
            break;
         }
         count++;
      }
    //  this->triangles->push_back(GbTriFaceMesh3D::TriFace(v1->id, v2->id, v3->id));
    // das geht bei Winkelplatte und Bathe
      this->triangles->push_back(GbTriFaceMesh3D::TriFace(id2, id1, id3));
    //  this->triangles->push_back(GbTriFaceMesh3D::TriFace(id1, id2, id3));
      countTris++;
   }

   std::cout<<"#################################"<<std::endl;
   std::cout<<"countTris:"<<countTris<<std::endl;
   std::cout<<"vecSize:"<<this->triangles->size()<<std::endl;
   this->attributes->resize(nodes->size());

   countTris=0;
   for(int u=0;u<(int)this->triangles->size(); u++)
   {
       double area = (*this->triangles)[u].getArea(*this->nodes);
       if(UbMath::zero(area))  countTris++;
   }
   std::cout<<"#################################"<<std::endl;
   std::cout<<"Area 0 für:"<<countTris<<" Dreiecke"<<std::endl;

   this->calculateValues();

   this->createVertexTriFaceMap();
}
/*===============================================================================*/
void FeAdhocTriFaceMesh3D::writeAdhoCMeshForStefan(string filename)
{
   std::cout << "FeAdhocTriFaceMesh3D::writeAdhoCMeshForStefan ...\n";
   List *elementlist = mesh->TL;

   ListItem* LI;
   Tri*      triangle;

   vector<GbPoint3D*>*     tmnodes     = new vector<GbPoint3D*>;
   vector<GbTriangle3D*>*  tmtriangles = new vector<GbTriangle3D*>;

   FOR_ALL(elementlist->first, LI, 0)
   {
      triangle = get_Tri(LI);

      GbPoint3D *node1 = new GbPoint3D(triangle->V[0]->x, triangle->V[0]->y, triangle->V[0]->z);
      GbPoint3D *node2 = new GbPoint3D(triangle->V[1]->x, triangle->V[1]->y, triangle->V[1]->z);
      GbPoint3D *node3 = new GbPoint3D(triangle->V[2]->x, triangle->V[2]->y, triangle->V[2]->z);

      tmnodes->push_back(node1);
      tmnodes->push_back(node2);
      tmnodes->push_back(node3);
      tmtriangles->push_back(new GbTriangle3D(node1, node2, node3));

   }

   GbTriangularMesh3D tmmesh("Name", tmnodes, tmtriangles);
   UbFileOutputASCII out(filename);
   tmmesh.writeAVSMesh(&out);
}

