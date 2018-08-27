#include "./QGbSphere3DInstrument.h"

/**** Qt ****/
#include <qlineedit.h>
#include <qstring.h>
#include <qcheckbox.h>

/**** CAB ****/
#include "./../GbSphere3D.h"

QGbSphere3DInstrument::QGbSphere3DInstrument( QWidget* parent, Qt::WFlags flags ):QDialog(parent,flags)
{

	ui.setupUi(this);
   
	this->gbSphere = NULL;
}

QGbSphere3DInstrument::~QGbSphere3DInstrument(void)
{
}

void QGbSphere3DInstrument::setGbSphere3D(GbSphere3D* sphere)
{
	this->gbSphere = sphere;
	ui.lineEditX->setText( QString("%1").arg(gbSphere->getX1Centroid() ) );
	ui.lineEditY->setText( QString("%1").arg(gbSphere->getX2Centroid() ) );
	ui.lineEditZ->setText( QString("%1").arg(gbSphere->getX3Centroid() ) );
   ui.lineEditName->setText( QString(gbSphere->getName().c_str()) );
	ui.lineEditRadius->setText( QString("%1").arg(gbSphere->getRadius() ) );
	ui.checkBoxActive->setChecked( true );
}

GbSphere3D* QGbSphere3DInstrument::getGbSphere3D(void)
{
	return this->gbSphere;
}

//void QGbSphere3DInstrument::SetGbObject3D(GbObject3D* gbObj)
//{
//		this->SetGbSphere(dynamic_cast<GbSphere3D*>(gbObj));
//}

void QGbSphere3DInstrument::on_pBtnOK_clicked()
{
	this->gbSphere->setCenterX1Coordinate(ui.lineEditX->text().toDouble());
	this->gbSphere->setCenterX2Coordinate(ui.lineEditY->text().toDouble());
	this->gbSphere->setCenterX3Coordinate(ui.lineEditZ->text().toDouble());
	this->gbSphere->setRadius(ui.lineEditRadius->text().toDouble());
   this->gbSphere->setName(ui.lineEditName->text().toStdString());
	//this->gbSphere->setActive( this->checkBoxActive->isChecked() );
	this->gbSphere->notifyObserversObjectChanged();
	this->accept();
}


void QGbSphere3DInstrument::on_pBtnCancel_clicked()
{
	this->reject();
}
