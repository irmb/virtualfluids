#include "./vtkGbCuboid3D.h"

#include "./../GbCuboid3D.h"
#include "./../../../userinterface/presentation/vtkEventCallbacks.h"

#include "vtkCubeSource.h"
#include "vtkPolyDataMapper.h"
//#include "math.h"

vtkGbCuboid3D::vtkGbCuboid3D(GbCuboid3D* gbObject)
{
	this->gbCuboid = gbObject;
	this->gbCuboid->addObserver(this);

   this->setName("vtkGbCuboid3D");

	this->source = vtkCubeSource::New();
   this->mapper = vtkPolyDataMapper::New();

	this->setValues();

	this->mapper->SetInput( this->source->GetOutput() );
	this->actor->SetMapper( this->mapper );
}

vtkGbCuboid3D::~vtkGbCuboid3D(void)
{
   this->gbCuboid->removeObserver(this);
	if (this->source) this->source->Delete();
}

//void vtkGbCuboid3D::ModifiedEventFired()
//{
//	//double a_orien[3];
//	double a_pos[3];
//	this->actor->GetPosition(a_pos);
//	//this->actor->GetOrientation(a_orien);
//	this->actor->SetPosition(0.0,0.0,0.0);
//	this->actor->SetOrientation(0.0,0.0,0.0);
//	this->actor->SetScale(1.0,1.0,1.0);
//
//	//cout<<"Orien:"<<a_orien[0]<<","<<a_orien[1]<<","<<a_orien[3]<<endl;
//	//cout<<"Position:"<<a_pos[0]<<","<<a_pos[1]<<","<<a_pos[3]<<endl;
//
//	this->gbCuboid->translate(a_pos[0], a_pos[1], a_pos[2]);
//	this->gbCuboid->notifyObserversObjectChanged();
//}

void vtkGbCuboid3D::applyActorModifications()
{
	if (isModified) 
	{
		double pos[3];
		double scale[3];
		//double orien[3];
		this->actor->GetPosition(pos);
		this->actor->GetScale(scale);
		//this->actor->GetOrientation(orien);

		this->actor->SetPosition(0.0,0.0,0.0);
		this->actor->SetOrientation(0.0,0.0,0.0);
		this->actor->SetScale(1.0,1.0,1.0);

		//cout<<"Orien:"<<a_orien[0]<<","<<a_orien[1]<<","<<a_orien[3]<<endl;
		//cout<<"Position:"<<a_pos[0]<<","<<a_pos[1]<<","<<a_pos[3]<<endl;


		////////////////////////////////////////////////////////////////////////////
		////Rotieren
		////[Cy x1 + Sy x3, x2, -Sy x1 + Cy x3, 1]
		//double center[3];
		//center[0] = this->gbCuboid->getX1Centroid();
		//center[1] = this->gbCuboid->getX2Centroid();
		//center[2] = this->gbCuboid->getX3Centroid();

		////Punkt1
		//double p1x = this->gbCuboid->getPoint1()->getX1Coordinate();
		//double p1y = this->gbCuboid->getPoint1()->getX2Coordinate();
		//double p1z = this->gbCuboid->getPoint1()->getX3Coordinate();

		//p1x = cos(orien[1]) * p1x + sin(orien[1]) * p1z;
		////p1y = p1y;
		//p1z = -sin(orien[1]) * p1x + cos(orien[1]) * p1z;

		//this->gbCuboid->getPoint1()->setX1(p1x);
		//this->gbCuboid->getPoint1()->setX2(p1y);
		//this->gbCuboid->getPoint1()->setX3(p1z);

		//
		////Punkt2
		//double p2x = this->gbCuboid->getPoint2()->getX1Coordinate();
		//double p2y = this->gbCuboid->getPoint2()->getX2Coordinate();
		//double p2z = this->gbCuboid->getPoint2()->getX3Coordinate();

		//p2x = cos(orien[1]) * p2x + sin(orien[1]) * p2z;
		////p1y = p1y;
		//p2z = -sin(orien[1]) * p2x + cos(orien[1]) * p2z;

		//this->gbCuboid->getPoint2()->setX1(p2x);
		//this->gbCuboid->getPoint2()->setX2(p2y);
		//this->gbCuboid->getPoint2()->setX3(p2z);
		//
		////////////////////////////////////////////////////////////////////////////

		if (scale[0] != 1.0) this->gbCuboid->scale(scale[0], scale[1], scale[2]);
		else this->gbCuboid->translate(pos[0], pos[1], pos[2]);
		this->gbCuboid->notifyObserversObjectChanged();

		//Methode der Basisklasse aufrufen.
		vtkPoElement3D::applyActorModifications();
	}
}

void vtkGbCuboid3D::setValues(void)
{
	double bounds[6];
	bounds[0] = this->gbCuboid->getX1Minimum();
	bounds[1] = this->gbCuboid->getX1Maximum();
	bounds[2] = this->gbCuboid->getX2Minimum();
	bounds[3] = this->gbCuboid->getX2Maximum();
	bounds[4] = this->gbCuboid->getX3Minimum();
	bounds[5] = this->gbCuboid->getX3Maximum();
	this->source->SetBounds(bounds);

//	this->actor->SetVisibility( this->gbCuboid->isActive() );
}

bool vtkGbCuboid3D::isPointInObject(double const point[3])
{
	return this->gbCuboid->isPointInGbObject3D(point[0], point[1], point[2]);
}

void vtkGbCuboid3D::objectChanged(UbObservable*)
{
	this->setValues();
	this->source->Update();
}

void vtkGbCuboid3D::objectWillBeDeleted(UbObservable*)
{
	//TODO: Hier muss auf jeden Fall noch was geschehen....
	this->gbCuboid->removeObserver(this);
	delete this;
}
