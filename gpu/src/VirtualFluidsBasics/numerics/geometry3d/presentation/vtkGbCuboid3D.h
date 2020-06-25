#ifndef VTKGBCUBOID3D_H
#define VTKGBCUBOID3D_H

#include "./../../../userinterface/presentation/vtkPoElement3D.h"

/**** vtk ****/
class vtkCubeSource;
class vtkPolyDataMapper;

class GbCuboid3D;

class vtkGbCuboid3D : public vtkPoElement3D
{
public:
	vtkGbCuboid3D(GbCuboid3D*);
	~vtkGbCuboid3D(void);
	void objectChanged(UbObservable*);
	void objectWillBeDeleted(UbObservable*);
	//void ModifiedEventFired(void);
	void applyActorModifications();
	bool isPointInObject(double const point[3]);
protected:
	void setValues();

	GbCuboid3D* gbCuboid;
	vtkCubeSource* source;
   vtkPolyDataMapper* mapper;
};
#endif   
