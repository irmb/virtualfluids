#ifndef VTKGBCYLINDER3D_H
#define VTKGBCYLINDER3D_H

#include "./../../../userinterface/presentation/vtkPoElement3D.h"

class GbCylinder3D;

class vtkCylinderSource;
class vtkPolyDataMapper;


class vtkGbCylinder3D : public vtkPoElement3D
{
public:
	vtkGbCylinder3D(GbCylinder3D* cylinder);
	~vtkGbCylinder3D();
	void objectChanged(UbObservable*);
	void objectWillBeDeleted(UbObservable*);
	//void ModifiedEventFired(void);
	void applyActorModifications();             
	bool isPointInObject(double const point[3]);
protected:
	void setValues();

	GbCylinder3D* gbCylinder;
   vtkCylinderSource* source;
   vtkPolyDataMapper* mapper;
};
#endif   

