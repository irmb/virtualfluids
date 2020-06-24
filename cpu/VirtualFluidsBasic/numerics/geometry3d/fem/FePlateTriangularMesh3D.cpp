#include <numerics/geometry3d/fem/FePlateTriangularMesh3D.h>

using namespace std;
                         
FePlateTriangularMesh3D::FePlateTriangularMesh3D() : GbTriangularMesh3D()
{
}

FePlateTriangularMesh3D::FePlateTriangularMesh3D(string name, vector<GbPoint3D*> *nodes, vector<GbTriangle3D*> *triangles) : GbTriangularMesh3D(name, nodes, triangles)
{
}

/*=============================================================================================*/

FePlateTriangularMesh3D::FePlateTriangularMesh3D(string name, vector<GbTriangle3D*> *triangles) : GbTriangularMesh3D(name, triangles)
{
}

/*=============================================================================================*/
FePlateTriangularMesh3D::FePlateTriangularMesh3D(string name, vector<GbPoint3D*> *nodes, vector<GbLine3D*> *edges, vector<GbTriangle3D*> *triangles) : GbTriangularMesh3D(name, nodes, edges, triangles)
{
}
/*=============================================================================================*/
FePlateTriangularMesh3D::~FePlateTriangularMesh3D()
{

}
/*======================================================================*/

/*======================================================================*/
bool FePlateTriangularMesh3D::isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
{
	//cout<<"GbTriangularMesh3D::isPointInGbObject3D"<<endl;
	//Sebastian
	//double xmin=this->getX1Minimum();	double xmax=this->getX1Maximum();
 //  double ymin=this->getX2Minimum();	double ymax=this->getX2Maximum();
 //  double zmin=this->getX3Minimum();	double zmax=this->getX3Maximum();
 //  double dX = (xmax-xmin)/100.;
 //  double dY = (ymax-ymin)/100.;
 //  double dZ = (zmax-zmin)/100.;
 //  GbCuboid3D boundingCube(xmin-dX, ymin-dY, zmin-dZ, xmax+dX, ymax+dY, zmax+dZ);
	//
	//if(this->isPointInObject3DHalfSpace(x1,x2,x3))
	//{
	//	return true;
	//}
	//if(!boundingCube.isPointInGbObject3D(x1, x2, x3))
	//{
	//	return false;
	//}

	//return false;
	//Marco
	int inside = 0;
   int nx1 = this->getNodesX1Dir()-1;
   int nx2 = this->getNodesX2Dir()-1;
   int maxTriangels = 2*nx1*nx2;

   //Überprüft, ob der Punkt innerhalb des Netzes liegt
	double xmin=this->getX1Minimum();	double xmax=this->getX1Maximum();
	double ymin=this->getX2Minimum();	double ymax=this->getX2Maximum();
	double zmin=this->getX3Minimum();	double zmax=this->getX3Maximum();
	if(	x1<=xmax && x1>=xmin
		&& x2<=ymax && x2>=ymin
		&& x3<=zmax && x3>=zmin)
	{
      //Achtung Sonderfall
      //Es wird nur gegen die obere Lage Dreiecke getestet, da die untere Lage um den Abstand 'height' verschoben ist!
      //Die Seiten müssen somit auch nicht berücksichtigt werden
      for(int i=0; i<int(maxTriangels);i++)     
		{
			if(	(triangles->at(i))->enclosesPoint2D(x1, x2)
				&&	x3<triangles->at(i)->getX3Centroid()
				&& x3>triangles->at(i)->getX3Centroid() - this->height  )
			{
					inside++;
			}
		}
	}

	if(inside!=0)
	{
		return true;
	}
	else return false;
}
/*======================================================================*/
FePlateTriangularMesh3D* FePlateTriangularMesh3D::createMeshByElements(int ElementsX1, int ElementsX2, double nulllevel, double deltaX1, double deltaX2, double height)
{
   FePlateTriangularMesh3D* triMesh = FePlateTriangularMesh3D::createMeshByNodes(ElementsX1+1, ElementsX2+1, nulllevel, deltaX1, deltaX2, height);
   return triMesh;
}
/*======================================================================*/
FePlateTriangularMesh3D* FePlateTriangularMesh3D::createMeshByNodes(int nodesX1, int nodesX2, double nulllevel, double deltaX1, double deltaX2, double height)
{
   cout<<"erstelle GbTriangularMesh3D -> ";
   vector<GbTriangle3D*>  *triangles = new vector<GbTriangle3D*>;
   vector<GbPoint3D*> *points = new vector<GbPoint3D*>;
   double h1=nulllevel+0.5*height;
   double h2=nulllevel-0.5*height;
   GbPoint3D* point1;
   GbPoint3D* point2;
   GbPoint3D* point3;
   GbPoint3D* point4;
   //Erstelle Knoten
   //oben
   for(int i=0; i<nodesX1; i++)
   {
      for(int j=0; j<nodesX2; j++)
      {
         GbPoint3D* point = new GbPoint3D(i*deltaX1,j*deltaX2,h1);
         points->push_back(point);
      }
   }
   //unten
   for(int i=0; i<nodesX1; i++)
   {
      for(int j=0; j<nodesX2; j++)
      {
         GbPoint3D* point = new GbPoint3D(i*deltaX1,j*deltaX2,h2);
         points->push_back(point);
      }
   }
   int size=int(points->size());
   //Erstelle Dreiecke
   //oben
   for(int i=0; i<(size/2)-nodesX2; i+=nodesX2)
   {
      for(int j=0; j<nodesX2-1; j++)
      {
         point1 = points->at(i+j);
         point2 = points->at(i+j+1);
         point3 = points->at(i+j+nodesX2);
         point4 = points->at(i+j+nodesX2+1);
         GbTriangle3D* tri1 = new GbTriangle3D(point1,point3,point4);
         GbTriangle3D* tri2 = new GbTriangle3D(point1,point4,point2);
         triangles->push_back(tri1);
         triangles->push_back(tri2);
      }
   }
   //unten
   for(int i=(size/2); i<size-nodesX2; i+=nodesX2)
   {
      for(int j=0; j<nodesX2-1; j++)
      {
         point1 = points->at(i+j);
         point2 = points->at(i+j+1);
         point3 = points->at(i+j+nodesX2);
         point4 = points->at(i+j+nodesX2+1);
         GbTriangle3D* tri1 = new GbTriangle3D(point4,point3,point1);
         GbTriangle3D* tri2 = new GbTriangle3D(point2,point4,point1);
         triangles->push_back(tri1);
         triangles->push_back(tri2);
      }
   }
   //Rand
   //Nord
   for(int i=0;i<nodesX1-1;i++)
   {	
      int a = i*nodesX2+nodesX2-1;
      int b = (i+1)*nodesX2+nodesX2-1;
      point1 = points->at(a);
      point2 = points->at(b);
      point3 = points->at(a+size/2);
      point4 = points->at(b+size/2);
      GbTriangle3D* tri1 = new GbTriangle3D(point1,point2,point3);
      GbTriangle3D* tri2 = new GbTriangle3D(point2,point4,point3);
      triangles->push_back(tri1);
      triangles->push_back(tri2);
   }
   //Süd
   for(int i=0;i<nodesX1-1;i++)
   {	
      int a = i*nodesX2;
      int b = (i+1)*nodesX2;
      point1 = points->at(a);
      point2 = points->at(b);
      point3 = points->at(a+size/2);
      point4 = points->at(b+size/2);
      GbTriangle3D* tri1 = new GbTriangle3D(point3,point2,point1);
      GbTriangle3D* tri2 = new GbTriangle3D(point3,point4,point2);
      triangles->push_back(tri1);
      triangles->push_back(tri2);
   }
   //Ost
   for(int j=0;j<nodesX2-1;j++)
   {	
      int a = j;
      int b = j+1;
      point1 = points->at(a);
      point2 = points->at(b);
      point3 = points->at(a+size/2);
      point4 = points->at(b+size/2);
      GbTriangle3D* tri1 = new GbTriangle3D(point1,point2,point3);
      GbTriangle3D* tri2 = new GbTriangle3D(point4,point3,point2);
      triangles->push_back(tri1);
      triangles->push_back(tri2);
   }
   //West
   for(int j=0;j<nodesX2-1;j++)
   {	
      int a = j+(nodesX1-1)*nodesX2;
      int b = j+(nodesX1-1)*nodesX2+1;
      point1 = points->at(a);
      point2 = points->at(b);
      point3 = points->at(a+size/2);
      point4 = points->at(b+size/2);
      GbTriangle3D* tri1 = new GbTriangle3D(point3,point2,point1);
      GbTriangle3D* tri2 = new GbTriangle3D(point2,point3,point4);
      triangles->push_back(tri1);
      triangles->push_back(tri2);
   }
   string name = "TriMesh";
   FePlateTriangularMesh3D* triMesh = new FePlateTriangularMesh3D(name,points,triangles);
   triMesh->setNullLevel((float)nulllevel);
   triMesh->setHeight(height);
   triMesh->setNodesX1Dir(nodesX1);
   triMesh->setNodesX2Dir(nodesX2);
   cout<<"done"<<endl;
   return triMesh;
}


