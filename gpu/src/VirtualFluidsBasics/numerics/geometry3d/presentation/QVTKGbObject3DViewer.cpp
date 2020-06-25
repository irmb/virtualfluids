#include "./QVTKGbObject3DViewer.h"

/**** Qt ****/
#include <qtabwidget.h>
#include <qlabel.h>
/**** vtk ****/

#include <QVTKWidget.h>
//#include "QvtkWindow.h"

/**** CAB ****/
#include "./../../../basics/utilities/UbMath.h"
#include "./../GbObject3DManager.h"


#include "./../../../userinterface/instrument/QManagerPresentatorInstrument.h"
//#include "./QGbObject3DManagerInstrument.h"



QVTKGbObject3DViewer::QVTKGbObject3DViewer():QVTKViewer3DApplication()
{
	//GbObjectManagerInstrument
   this->gbObject3DManager = new GbObject3DManager();
   QManagerPresentatorInstrument* gbObjManInst = new QManagerPresentatorInstrument(gbObject3DManager);
	//gbObjManInst->setQViewer(this->getViewer());
	
	//Instrumente hinzufügen
	this->addInstrument(gbObjManInst, "Geometries");
}

QVTKGbObject3DViewer::~QVTKGbObject3DViewer()
{
}

