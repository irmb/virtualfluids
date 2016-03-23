#ifndef QVTKGBOBJECT3DVIEWER_H
#define QVTGBOBJECT3DKVIEWER_H

#include "./../../../userinterface/viewer3d/QVTKViewer3DApplication.h"

class QVTKWindow;
class QVTKViewer3D;
class GbObject3DManager;
class OctNodeGridManager;

class QVTKGbObject3DViewer : public QVTKViewer3DApplication
{
public:
	QVTKGbObject3DViewer();
	~QVTKGbObject3DViewer();
   
protected:

	GbObject3DManager* gbObject3DManager;

};
#endif
