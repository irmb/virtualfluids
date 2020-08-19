#include "./QGbObject3DInstrument.h"

/**** Qt ****/
#include <qlineedit.h>

/**** vtk ****/
#include "./../GbObject3D.h"
#include "./../../../basics/utilities/UbMath.h"

//#define PI   3.14159265358979323846

QGbObject3DInstrument::QGbObject3DInstrument( QWidget* parent, Qt::WFlags flags )
{
	ui.setupUi(this);

	this->gbObject3D = NULL;
}

QGbObject3DInstrument::~QGbObject3DInstrument()
{
}

void QGbObject3DInstrument::setGbObject3D(GbObject3D* obj)
{                               
	this->gbObject3D = obj;
}

GbObject3D* QGbObject3DInstrument::getGbObject3D()
{
	return this->gbObject3D;
}

void QGbObject3DInstrument::on_pBtnOK_clicked()
{
	double rx = ui.lineEditRotationX->text().toDouble();
	double ry = ui.lineEditRotationY->text().toDouble();
	double rz = ui.lineEditRotationZ->text().toDouble();

	rx *= UbMath::PI /180;     
	ry *= UbMath::PI /180;
	rz *= UbMath::PI /180;

	if ( rx != 0.0 || ry != 0.0 || rz != 0.0 ) this->gbObject3D->rotate(rx, ry, rz);

	double sx = ui.lineEditScalingX->text().toDouble();
	double sy = ui.lineEditScalingY->text().toDouble();
	double sz = ui.lineEditScalingZ->text().toDouble();

	if ( sx != 0.0 || sy != 0.0 || sz != 0.0 ) this->gbObject3D->scale(sx, sy, sz);

	double x = ui.lineEditTranlationX->text().toDouble();
	double y = ui.lineEditTranlationY->text().toDouble();
	double z = ui.lineEditTranlationZ->text().toDouble();

	if ( x != 0.0 || y != 0.0 || z != 0.0 ) this->gbObject3D->translate(x, y, z);

	this->gbObject3D->notifyObserversObjectChanged();

	this->accept();
}


void QGbObject3DInstrument::on_pBtnCancel_clicked()
{
	this->reject();
}
