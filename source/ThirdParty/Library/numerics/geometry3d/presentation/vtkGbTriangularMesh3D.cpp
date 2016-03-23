#include "./vtkGbTriangularMesh3D.h"

/**** CAB ****/
#include "./../GbTriangularMesh3D.h"
#include "./../GbTriangle3D.h"
#include "./../../../basics/utilities/UbMath.h"

/**** vtk ****/
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkIdList.h"

/*** temp ****/
#include "./../GbPoint3D.h"
#include "vtkMatrix4x4.h"

//#define PI   3.14159265358979323846

vtkGbTriangularMesh3D::vtkGbTriangularMesh3D(GbTriangularMesh3D* mesh)
{
	this->gbTriangularMesh = mesh;
	this->gbTriangularMesh->addObserver(this);

   this->setName("vtkGbTriangularMesh3D");

	this->unstGrid = vtkUnstructuredGrid::New();
	this->buildGrid();

	
	this->dataSetMapper = vtkDataSetMapper::New();
	this->dataSetMapper->SetInput(unstGrid);

	this->actor->SetMapper( this->dataSetMapper );
}

vtkGbTriangularMesh3D::~vtkGbTriangularMesh3D(void)
{
   this->gbTriangularMesh->removeObserver(this);
	this->unstGrid->Delete();
	this->dataSetMapper->Delete();
}

//void vtkGbTriangularMesh3D::ModifiedEventFired()
//{
//	double pos[3];
//	this->actor->GetPosition(pos);
//
//	this->actor->SetPosition(0.0,0.0,0.0);
//	this->actor->SetOrientation(0.0,0.0,0.0);
//	this->actor->SetScale(1.0,1.0,1.0);
//
//	double x1 = pos[0];
//	double x2 = pos[1];
//	double x3 = pos[3];
//
//	vector<GbPoint3D*>* pointList = this->gbTriangularMesh->getNodes();
//	for (int pos=0; pos<pointList->size(); pos++) 
//	{
//		(*pointList)[pos]->translate(x1,x2,x3);
//		//((*pointList)[pos])->translate(pos[0], pos[1], pos[3]);
//	}
//	this->gbTriangularMesh->notifyObserversObjectChanged();
//}

void vtkGbTriangularMesh3D::applyActorModifications()
{
	if (isModified) 
	{
		double pos[3];
		double orien[3];
		double scale[3];
		this->actor->GetPosition(pos);
		this->actor->GetOrientation(orien);
		this->actor->GetScale(scale);

      orien[0] = orien[0] / 180 * UbMath::PI;
		orien[1] = orien[1] / 180 * UbMath::PI;
		orien[2] = orien[2] / 180 * UbMath::PI;

		//cout<<"Orien:"<<orien[0]<<","<<orien[1]<<","<<orien[3]<<endl;
		//cout<<"Position:"<<pos[0]<<","<<pos[1]<<","<<pos[3]<<endl;
		//cout<<"Scale:"<<scale[0]<<","<<scale[1]<<","<<scale[3]<<endl;
		
		this->actor->SetPosition(0.0,0.0,0.0);
		this->actor->SetOrientation(0.0,0.0,0.0);
		this->actor->SetScale(1.0,1.0,1.0);
		
		vector<GbPoint3D*>* pointList = this->gbTriangularMesh->getNodes();
		for (int index=0; index<(int)pointList->size(); index++) 
		{
			(*pointList)[index]->rotate(orien[0], orien[1], orien[2]);
			(*pointList)[index]->scale(scale[0], scale[1], scale[2]);
			(*pointList)[index]->translate(pos[0], pos[1], pos[2]);
		}
		this->gbTriangularMesh->notifyObserversObjectChanged();
		//Methode der Basisklasse aufrufen.
		vtkPoElement3D::applyActorModifications();
	}
}

void vtkGbTriangularMesh3D::buildGrid(void)
{
	this->unstGrid->Reset();

	vector<GbTriangle3D*>* triangles = this->gbTriangularMesh->getTriangles();
	double xyz[3];
	//this.setContext(new PoContext3D());

	vtkPoints* points  = vtkPoints::New();
	vtkTriangle* triangle = vtkTriangle::New();
	for(int u=0; u<(int)triangles->size(); u++)
	{
		xyz[0] = (*triangles)[u]->getPoint(0)->getX1Coordinate();
		xyz[1] = (*triangles)[u]->getPoint(0)->getX2Coordinate();
		xyz[2] = (*triangles)[u]->getPoint(0)->getX3Coordinate();
		triangle->GetPointIds()->InsertId(0, points->InsertNextPoint(xyz));
		//points.InsertPoint(u, xyz);       // 3D geometry

		xyz[0] = (*triangles)[u]->getPoint(1)->getX1Coordinate();
		xyz[1] = (*triangles)[u]->getPoint(1)->getX2Coordinate();
		xyz[2] = (*triangles)[u]->getPoint(1)->getX3Coordinate();
		triangle->GetPointIds()->InsertId(1, points->InsertNextPoint(xyz));
		//points.InsertPoint(u, xyz);       // 3D geometry

		xyz[0] = (*triangles)[u]->getPoint(2)->getX1Coordinate();
		xyz[1] = (*triangles)[u]->getPoint(2)->getX2Coordinate();
		xyz[2] = (*triangles)[u]->getPoint(2)->getX3Coordinate();
		triangle->GetPointIds()->InsertId(2, points->InsertNextPoint(xyz));
		//points.InsertPoint(u, xyz);       // 3D geometry

		this->insertNextCell( triangle->GetCellType(), triangle->GetPointIds() ); // grid topology

	}
	this->setPoints(points);
	//this->source->SetCenter(	this->gbSphere->getX1Centroid(),
	//	this->gbSphere->getX2Centroid(),
	//	this->gbSphere->getX3Centroid()	);

	//this->source->SetRadius( this->gbSphere->getRadius() );
	//this->actor->SetVisibility( this->gbSphere->isActive() );
	//this->unstGrid->Modified();
}

int vtkGbTriangularMesh3D::insertNextCell(int type, vtkIdList* idList)
{
	return this->unstGrid->InsertNextCell(type, idList);
}

void vtkGbTriangularMesh3D::setPoints(vtkPoints* points)
{
	this->unstGrid->SetPoints(points);
}

void vtkGbTriangularMesh3D::objectChanged(UbObservable*)
{
	this->buildGrid();
	this->unstGrid->Update();
}

void vtkGbTriangularMesh3D::objectWillBeDeleted(UbObservable*)
{
	//TODO: Hier muss auf jeden Fall noch was geschehen....
	this->gbTriangularMesh->removeObserver(this);
	delete this;
}
