#include "./QGbCylinder3DInstrument.h"

/**** Qt ****/
#include <QtCore/QString>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QTabWidget>

/**** CAB ****/
#include "./../GbCylinder3D.h"
#include "./../GbPoint3D.h"

QGbCylinder3DInstrument::QGbCylinder3DInstrument( QWidget* parent, Qt::WFlags flags )
{
	ui.setupUi(this);
   
   /*JZ TODO daher Buttons noch ausgeschaltet (29.11.05)*/
   ui.rBtnXAxis->setEnabled(false);
   ui.rBtnYAxis->setEnabled(false);
   ui.rBtnZAxis->setEnabled(false);
	
   this->gbCylinder = NULL;
}

QGbCylinder3DInstrument::~QGbCylinder3DInstrument(void)
{
}

void QGbCylinder3DInstrument::setGbCylinder3D(GbCylinder3D* cylinder)
{
   this->gbCylinder = cylinder;
   ui.lineEdit1_X_1->setText( QString("%1").arg(gbCylinder->getPoint1()->x1 ) );
   ui.lineEdit1_Y_1->setText( QString("%1").arg(gbCylinder->getPoint1()->x2 ) );
   ui.lineEdit1_Z_1->setText( QString("%1").arg(gbCylinder->getPoint1()->x3 ) );
   ui.lineEdit1_X_2->setText( QString("%1").arg(gbCylinder->getPoint2()->x1 ) );
   ui.lineEdit1_Y_2->setText( QString("%1").arg(gbCylinder->getPoint2()->x2 ) );
   ui.lineEdit1_Z_2->setText( QString("%1").arg(gbCylinder->getPoint2()->x3 ) );
   ui.dSpBoxRadius1->setValue(gbCylinder->getRadius());
   ui.checkBoxActive1->setChecked( true );
   ui.lineEdit2_X->setText( QString("%1").arg(gbCylinder->getPoint1()->x1 ) );
   ui.lineEdit2_Y->setText( QString("%1").arg(gbCylinder->getPoint1()->x2 ) );
   ui.lineEdit2_Z->setText( QString("%1").arg(gbCylinder->getPoint1()->x3 ) );
   ui.dSpBoxRadius2->setValue(gbCylinder->getRadius());
   ui.checkBoxActive2->setChecked( true );
   ui.lineEditLength->setText( QString("%1").arg(gbCylinder->getHeight()) );
   //if (!this->gbCylinder->isParallelToX1Axis()) 
   //{
   //   if (!this->gbCylinder->isParallelToX2Axis()) 
   //   {
   //      ui.rBtnZAxis->setChecked(true);
   //   }
   //   else ui.rBtnYAxis->setChecked(true);
   //}
   //else ui.rBtnXAxis->setChecked(true);
}

GbCylinder3D* QGbCylinder3DInstrument::getGbCylinder3D(void)
{
	return this->gbCylinder;
}

//void QGbSphere3DInstrument::SetGbObject3D(GbObject3D* gbObj)
//{
//		this->SetGbSphere(dynamic_cast<GbSphere3D*>(gbObj));
//}

void QGbCylinder3DInstrument::on_pBtnOK_clicked()
{
   if(ui.tabWidget->currentIndex()==0)
   {
      this->gbCylinder->setPoint1(  ui.lineEdit1_X_1->text().toDouble(),
                                    ui.lineEdit1_Y_1->text().toDouble(),
                                    ui.lineEdit1_Z_1->text().toDouble());
     
      this->gbCylinder->setPoint2(  ui.lineEdit1_X_2->text().toDouble(),
                                    ui.lineEdit1_Y_2->text().toDouble(),
                                    ui.lineEdit1_Z_2->text().toDouble());
      this->gbCylinder->setRadius(ui.dSpBoxRadius1->value());

      this->gbCylinder->notifyObserversObjectChanged();
   }
   if(ui.tabWidget->currentIndex()==1)
   {
      this->gbCylinder->setPoint1(  ui.lineEdit2_X->text().toDouble(),
                                    ui.lineEdit2_Y->text().toDouble(),
                                    ui.lineEdit2_Z->text().toDouble());
      this->gbCylinder->setPoint2(  ui.lineEdit2_X->text().toDouble(),
                                    ui.lineEdit2_Y->text().toDouble()+ui.lineEditLength->text().toDouble(),
                                    ui.lineEdit2_Z->text().toDouble());
      this->gbCylinder->setRadius(ui.dSpBoxRadius2->value());

      this->gbCylinder->notifyObserversObjectChanged();
   }

   this->accept();
}


void QGbCylinder3DInstrument::on_pBtnCancel_clicked()
{
	this->reject();
}
