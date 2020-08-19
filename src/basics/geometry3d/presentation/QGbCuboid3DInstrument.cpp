#include "./QGbCuboid3DInstrument.h"

/**** Qt ****/
#include <qlineedit.h>
#include <qstring.h>
#include <qcheckbox.h>

/**** CAB ****/
#include "./../GbCuboid3D.h"
#include "./../GbPoint3D.h"


QGbCuboid3DInstrument::QGbCuboid3DInstrument( QWidget* parent, Qt::WFlags flags ):QDialog(parent, flags)
{
	ui.setupUi(this);

	this->gbCuboid = NULL;

}

QGbCuboid3DInstrument::~QGbCuboid3DInstrument()
{
}

void QGbCuboid3DInstrument::setGbCuboid3D(GbCuboid3D* cuboid)
{
	this->gbCuboid = cuboid;
	ui.lineEditPoint1X->setText( QString("%1").arg(gbCuboid->getPoint1()->getX1Coordinate() ) );
	ui.lineEditPoint1Y->setText( QString("%1").arg(gbCuboid->getPoint1()->getX2Coordinate() ) );
	ui.lineEditPoint1Z->setText( QString("%1").arg(gbCuboid->getPoint1()->getX3Coordinate() ) );
	ui.lineEditPoint2X->setText( QString("%1").arg(gbCuboid->getPoint2()->getX1Coordinate() ) );
	ui.lineEditPoint2Y->setText( QString("%1").arg(gbCuboid->getPoint2()->getX2Coordinate() ) );
	ui.lineEditPoint2Z->setText( QString("%1").arg(gbCuboid->getPoint2()->getX3Coordinate() ) );
	//this->checkBoxActive->setChecked( cuboid->isActive() );
	ui.checkBoxActive->setChecked( true );
}

GbCuboid3D* QGbCuboid3DInstrument::getGbCuboid3D()
{
	return this->gbCuboid;
}

//void QGbCuboid3DInstrument::SetGbObject3D(GbObject3D* gbObj)
//{
//		this->SetGbSphere(dynamic_cast<GbSphere3D*>(gbObj));
//}

void QGbCuboid3DInstrument::on_pBtnOK_clicked()
{
	this->gbCuboid->getPoint1()->setX1(ui.lineEditPoint1X->text().toDouble() );
	this->gbCuboid->getPoint1()->setX2(ui.lineEditPoint1Y->text().toDouble() );
	this->gbCuboid->getPoint1()->setX3(ui.lineEditPoint1Z->text().toDouble() );
	this->gbCuboid->getPoint2()->setX1(ui.lineEditPoint2X->text().toDouble() );
	this->gbCuboid->getPoint2()->setX2(ui.lineEditPoint2Y->text().toDouble() );
	this->gbCuboid->getPoint2()->setX3(ui.lineEditPoint2Z->text().toDouble() );
	//this->gbCuboid->setActive( this->checkBoxActive->isChecked() );

	this->gbCuboid->notifyObserversObjectChanged();
	this->accept();
}


void QGbCuboid3DInstrument::on_pBtnCancel_clicked()
{
	this->reject();
}

