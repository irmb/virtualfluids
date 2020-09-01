#include "./vtkGbSphere3D.h"

#include "./../GbSphere3D.h"
#include "./../../../userinterface/presentation/vtkEventCallbacks.h"

#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"


vtkGbSphere3D::vtkGbSphere3D(GbSphere3D* gbObject):vtkPoElement3D()
{
	this->gbSphere = gbObject;
	this->gbSphere->addObserver(this);
   
   this->setName("vtkGbSphere3D");


	this->source = vtkSphereSource::New();
   this->mapper = vtkPolyDataMapper::New();
	
	this->setValues();

	this->mapper->SetInput( this->source->GetOutput() );
	this->actor->SetMapper( this->mapper );
//	this->actor->GetProperty()->SetRepresentationToWireframe();
}

vtkGbSphere3D::~vtkGbSphere3D(void)
{
   this->gbSphere->removeObserver(this);
	if (this->source) this->source->Delete();
}

//void vtkGbSphere3D::ModifiedEventFired()
//{
//	////double a_orien[3];
//	//double a_pos[3];
//	////double a_scale[3];
//	//this->actor->GetPosition(a_pos);
//	////this->actor->GetOrientation(a_orien);
//	////this->actor->GetScale(a_scale);
//
//	//this->actor->SetPosition(0.0,0.0,0.0);
//	//this->actor->SetOrientation(0.0,0.0,0.0);
//	//this->actor->SetScale(1.0,1.0,1.0);
//	//
//	////cout<<"Orien:"<<a_orien[0]<<","<<a_orien[1]<<","<<a_orien[3]<<endl;
//	////cout<<"Position:"<<a_pos[0]<<","<<a_pos[1]<<","<<a_pos[3]<<endl;
//	////cout<<"Scale:"<<a_scale[0]<<","<<a_scale[1]<<","<<a_scale[3]<<endl;
//
//	//this->gbSphere->translate(a_pos[0], a_pos[1], a_pos[2]);
//	//this->gbSphere->notifyObserversObjectChanged();
//	PoElement3D::ModifiedEventFired();
//}

void vtkGbSphere3D::setValues(void)
{
	this->source->SetCenter(	this->gbSphere->getX1Centroid(),
								this->gbSphere->getX2Centroid(),
								this->gbSphere->getX3Centroid()	);
	
	this->source->SetRadius( this->gbSphere->getRadius() );
//	this->actor->SetVisibility( this->gbSphere->isActive() );
}

void vtkGbSphere3D::applyActorModifications()
{
	if (this->isModified)
	{
		//double a_orien[3];
		double pos[3];
		double scale[3];
		this->actor->GetPosition(pos);
		//this->actor->GetOrientation(a_orien);
		this->actor->GetScale(scale);

		this->actor->SetPosition(0.0,0.0,0.0);
		this->actor->SetOrientation(0.0,0.0,0.0);
		this->actor->SetScale(1.0,1.0,1.0);

		//cout<<"Orien:"<<a_orien[0]<<","<<a_orien[1]<<","<<a_orien[3]<<endl;
		//cout<<"Position:"<<a_pos[0]<<","<<a_pos[1]<<","<<a_pos[3]<<endl;
		//cout<<"Scale:"<<a_scale[0]<<","<<a_scale[1]<<","<<a_scale[3]<<endl;

		if (scale[0] != 1.0) this->gbSphere->scale(scale[0], scale[1], scale[2]);
		else this->gbSphere->translate(pos[0], pos[1], pos[2]);
		this->gbSphere->notifyObserversObjectChanged();

		vtkPoElement3D::applyActorModifications();
	}
}

bool vtkGbSphere3D::isPointInObject(double const point[3])
{
	return this->gbSphere->isPointInGbObject3D(point[0], point[1], point[2]);
}

//Wird aufgerufen, wenn sich das zugehörige GBObject3D ändert.
void vtkGbSphere3D::objectChanged(UbObservable*)
{
	this->setValues();
	this->source->Modified();
	this->actor->Modified();
}

void vtkGbSphere3D::objectWillBeDeleted(UbObservable*)
{
	//TODO: Hier muss auf jeden Fall noch was geschehen....
	this->gbSphere->removeObserver(this);
	delete this;
}

