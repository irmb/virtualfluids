#ifndef STL2INP_H
#define STL2INP_H

#include <QtGui/QMainWindow>
#include <QtGui/QProgressBar>
#include <QTimer>
#include "stl2inpUI.h"

class STL2INP : public QMainWindow
{
	Q_OBJECT

public:
	STL2INP(QWidget *parent = 0, Qt::WFlags flags = 0);
	~STL2INP();
	
	int startLevel, stopLevel;
	double delta;
	int nx1, nx2, nx3;

private:
	Ui::STL2INPClass ui;

	private slots:
		void on_checkBox_stateChanged(int);
		void on_pBtn_Input_pressed();
		void on_pBtn_Output_pressed();
		void on_pBtn_Output_2_pressed();
		void on_pBtn_Convert_pressed();
		void on_lineEdit_In_returnPressed();
		void on_lineEdit_Out_returnPressed();
		void on_lineEdit_Out_2_returnPressed();
		void on_pBtn_EditMesh_pressed();
};

#endif // STL2INP_H
