#ifndef QGBCUBOID3DINSTRUMENT_H
#define QGBCUBOID3DINSTRUMENT_H


#include "./QGbCuboid3DInstrumentUI.h"

class GbCuboid3D;
class GbObject3D;

class QGbCuboid3DInstrument : public QDialog
{
	Q_OBJECT

public:
	QGbCuboid3DInstrument( QWidget* parent = 0, Qt::WFlags fl = 0 );
	~QGbCuboid3DInstrument();
	void setGbCuboid3D(GbCuboid3D* cuboid);        
	GbCuboid3D* getGbCuboid3D();

protected:
	GbCuboid3D* gbCuboid;

private:
	Ui::QGbCuboid3DInstrument ui;

private slots:
	void on_pBtnOK_clicked();
	void on_pBtnCancel_clicked();
};

#endif   
