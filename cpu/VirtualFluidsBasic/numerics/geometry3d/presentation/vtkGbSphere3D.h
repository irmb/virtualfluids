#ifndef VTKGBSPHERE3D_H
#define VTKGBSPHERE3D_H

#include "./../../../userinterface/presentation/vtkPoElement3D.h"
//#include "./../../../../vtkEventListeners.h"

class GbSphere3D;

class vtkSphereSource;
class vtkPolyDataMapper;

class vtkGbSphere3D : public vtkPoElement3D
{
public:
	vtkGbSphere3D(GbSphere3D*);
	~vtkGbSphere3D(void);
	void objectChanged(UbObservable*);
	void objectWillBeDeleted(UbObservable*);
	//void ModifiedEventFired(void);
	void applyActorModifications(); 
	bool isPointInObject(double const point[3]);

   virtual string toString() { return "vtkGbSphere3D";  }

protected:
	void setValues();

	GbSphere3D* gbSphere;
   vtkPolyDataMapper* mapper;
	vtkSphereSource* source;
};
#endif   

