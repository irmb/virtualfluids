#ifndef FEPLATETRIANGULARMESH_H
#define FEPLATETRIANGULARMESH_H

#include <sstream>
#include <iostream>

#include <geometry3d/GbTriangularMesh3D.h>


/*=========================================================================*/
/* GbTriangularMesh3D                                                      */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
*/
class FePlateTriangularMesh3D : public GbTriangularMesh3D 
{                             
public:
   FePlateTriangularMesh3D();
   FePlateTriangularMesh3D(std::string name, std::vector<GbPoint3D*> *nodes, std::vector<GbTriangle3D*> *triangles);
   FePlateTriangularMesh3D(std::string name, std::vector<GbTriangle3D*> *triangles);
   FePlateTriangularMesh3D(std::string name, std::vector<GbPoint3D*> *nodes, std::vector<GbLine3D*> *edges, std::vector<GbTriangle3D*> *triangles);
	virtual ~FePlateTriangularMesh3D();   

   static FePlateTriangularMesh3D* createMeshByNodes(int nodesX1, int nodesX2, double nulllevel, double deltaX1, double deltaX2, double height);
   static FePlateTriangularMesh3D* createMeshByElements(int ElementsX1, int ElementsX2, double nulllevel, double deltaX1, double deltaX2, double height);


	float	 getNullLevel()						{	return this->nulllevel;	}
	void	 setNullLevel(float nulllevel)	{	this->nulllevel = nulllevel;	}
   double getHeight()							{	return this->height;	}
   void	 setHeight(double height)			{	this->height = height;	}
   int    getNodesX1Dir()						{	return this->nodesX1Dir;	}
   void	 setNodesX1Dir(int nodesX1Dir)   {	this->nodesX1Dir = nodesX1Dir;	}
   int    getNodesX2Dir()						{	return this->nodesX2Dir;	}
   void	 setNodesX2Dir(int nodesX2Dir)   {	this->nodesX2Dir = nodesX2Dir;	}

	bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3);
   
   /*======================================================================*/
private:
	float			nulllevel;
	double		height;
   int         nodesX1Dir;
   int         nodesX2Dir;
};
/*=========================================================================*/
#endif
