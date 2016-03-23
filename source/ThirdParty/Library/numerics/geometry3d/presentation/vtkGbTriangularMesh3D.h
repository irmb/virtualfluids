#ifndef VTKGBTRIANGULARMESH3D_H
#define VTKGBTRIANGULARMESH3D_H

#include "./../../../userinterface/presentation/vtkPoElement3D.h"

class GbTriangularMesh3D;
class vtkUnstructuredGrid;
class vtkDataSetMapper;
class vtkIdList;
class vtkPoints;

class vtkGbTriangularMesh3D : public vtkPoElement3D
{
public:
	vtkGbTriangularMesh3D(GbTriangularMesh3D* mesh);
	~vtkGbTriangularMesh3D();
	void objectChanged(UbObservable* );
	void objectWillBeDeleted(UbObservable* );
	int insertNextCell(int, vtkIdList*);
	void setPoints(vtkPoints*);
	//void ModifiedEventFired(void);
	void applyActorModifications();
protected:
	void buildGrid();

	GbTriangularMesh3D* gbTriangularMesh;
	vtkUnstructuredGrid* unstGrid;
	vtkDataSetMapper* dataSetMapper;
};
#endif   

