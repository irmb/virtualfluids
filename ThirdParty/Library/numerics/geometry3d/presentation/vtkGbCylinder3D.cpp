#include "./vtkGbCylinder3D.h"

#include "./../GbCylinder3D.h"
#include "./../GbPoint3D.h"
#include "./../../../userinterface/presentation/vtkEventCallbacks.h"

#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"


vtkGbCylinder3D::vtkGbCylinder3D(GbCylinder3D* gbObject)
{
	this->gbCylinder = gbObject;
	this->gbCylinder->addObserver(this);

   this->setName("vtkGbCylinder3D");

	this->source = vtkCylinderSource::New();
   this->mapper = vtkPolyDataMapper::New();
	
	this->setValues();

	this->mapper->SetInput( this->source->GetOutput() );
	this->actor->SetMapper( this->mapper );

   //this->applyActorModifications();
}

vtkGbCylinder3D::~vtkGbCylinder3D(void)
{
   this->gbCylinder->removeObserver(this);
	if (this->source) this->source->Delete();
}


void vtkGbCylinder3D::setValues(void)
{
   //this->source->SetCenter(   this->gbCylinder->getX1Centroid(),
   //                           this->gbCylinder->getX2Centroid(),
   //                           this->gbCylinder->getX3Centroid());
   //this->source->SetHeight(this->gbCylinder->getLength());
   //this->source->SetRadius( this->gbCylinder->getRadius());

   /* JZ Attention not ready still some work TODO*/
   this->source->SetHeight(this->gbCylinder->getHeight());
   this->source->SetCenter(this->gbCylinder->getX1Centroid(),
                           this->gbCylinder->getX2Centroid(),
                           this->gbCylinder->getX3Centroid());
   this->source->SetRadius( this->gbCylinder->getRadius() );
   this->source->SetResolution(10);
}

void vtkGbCylinder3D::applyActorModifications()
{
   //this->actor->SetScale(1.0, this->gbCylinder->getLength(), 1.0);
   this->source->SetHeight(this->gbCylinder->getHeight());
   this->actor->SetPosition(  this->gbCylinder->getPoint1()->x1,
                              this->gbCylinder->getPoint1()->x2,
                              this->gbCylinder->getPoint1()->x3);
   this->source->SetRadius( this->gbCylinder->getRadius() );



   //if (this->isModified)
	//{
	//	double pos[3];
	//	double scale[3];
	//	this->actor->GetPosition(pos);
	//	this->actor->GetScale(scale);

	//	this->actor->SetPosition(0.0,0.0,0.0);
	//	this->actor->SetOrientation(0.0,0.0,0.0);
	//	this->actor->SetScale(1.0,1.0,1.0);


	//	if (scale[0] != 1.0) this->gbCylinder->scale(scale[0], scale[1], scale[2]);
	//	else this->gbCylinder->translate(pos[0], pos[1], pos[2]);
	//	this->gbCylinder->notifyObserversObjectChanged();

	//	vtkPoElement3D::applyActorModifications();
	//}
}

bool vtkGbCylinder3D::isPointInObject(double const point[3])
{
	return this->gbCylinder->isPointInGbObject3D(point[0], point[1], point[2]);
}

//Wird aufgerufen, wenn sich das zugehörige GBObject3D ändert.
void vtkGbCylinder3D::objectChanged(UbObservable*)
{
   this->setValues();
//	this->applyActorModifications();
	this->source->Modified();
}

void vtkGbCylinder3D::objectWillBeDeleted(UbObservable*)
{
	//TODO: Hier muss auf jeden Fall noch was geschehen....
	this->gbCylinder->removeObserver(this);
	delete this;
}

