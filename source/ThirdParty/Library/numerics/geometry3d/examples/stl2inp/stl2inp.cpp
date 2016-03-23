#include "stl2inp.h"
#include <QtGui/QFileDialog>
#include <QString>
#include <QFile>
#include <QMessageBox>
#include <cstdio>

#include "./QDefineUniformMesh.h"

#include "./../../../../../source/basics/utilities/UbFileInputASCII.h"
#include "./../../../../../source/basics/utilities/UbFileOutputASCII.h"
#include "./../../../../../source/basics/utilities/UbFileOutputBinary.h"
#include "./../../../../../source/numerics/geometry3d/GbTriangularMesh3D.h"
#include "./../../../../../source/numerics/geometry3d/creator/GbTriangularMesh3DCreator.h"
#include "./../../../../../source/numerics/geometry3D/CoordinateTransformation3D.h"
#include "./../../../../../source/basics/utilities/UbTiming.h"
#include "./../../../../../source/octree/facette/OctFacettenGrid2.h"

STL2INP::STL2INP(QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags)
{
	ui.setupUi(this);
	startLevel	= 0;
	stopLevel	= 3;
	delta		= 10.00;
	nx1			= 30;
	nx2			= 15;
	nx3			= 5;
	
}

STL2INP::~STL2INP()
{

}

void STL2INP::on_pBtn_Input_pressed()
{
	QString s = QFileDialog::getOpenFileName(
		this,
		"Choose a file",
		"E:/",
		"STL-Files (*.stl)");
	if(s != ""){
		ui.lineEdit_In->setText(s);
		ui.statusBar->showMessage("Input-File: Filename defined", 3000);
	}
	else
		ui.statusBar->showMessage("Input-File: No file found", 3000);

}

void STL2INP::on_lineEdit_In_returnPressed(){
	QString s = ui.lineEdit_In->text();
	if(s != ""){
		if(!s.endsWith(".stl",Qt::CaseSensitivity(false)))
		{
			s.append(".stl");
			ui.lineEdit_In->setText(s);
		}
		if(QFile::exists(s))
			ui.statusBar->showMessage("Inputput-File: File found", 3000);
		else
			ui.statusBar->showMessage("Input-File: File does not exist", 3000);
	}
	else
		ui.statusBar->showMessage("Input-File: no Filename", 3000);
}

void STL2INP::on_pBtn_Output_pressed(){
	QString s = QFileDialog::getSaveFileName(
		this,
		"Choose a filename to save under",
		"E:/",
		"AVS-File (*.inp)");
	if(s != ""){
		ui.lineEdit_Out->setText(s);
		ui.statusBar->showMessage("Output-File: Filename defined", 3000);
	}
	else
		ui.statusBar->showMessage("Output-File: No file defined", 3000);
}

void STL2INP::on_lineEdit_Out_returnPressed(){
	QString s = ui.lineEdit_Out->text();
	if(s != ""){
		if(!s.endsWith(".inp",Qt::CaseSensitivity(false)))
		{
			s.append(".inp");
			ui.lineEdit_Out->setText(s);
		}
		if (QFile::exists(s))
			if(QMessageBox::question(this,
				tr("Overwrite File? -- Application Name"),
				tr("A file called %1 already exists."
				"Do you want to overwrite it?")
				.arg(s),
				tr("&Yes"), tr("&No"),
				QString(), 0, 1))
				ui.lineEdit_Out->setText("");
			else
				ui.statusBar->showMessage("Output-File: overwrite existing File", 3000);
		else
			ui.statusBar->showMessage("Output-File: Filename defined", 3000);
	}
	else
		ui.statusBar->showMessage("Output-File: No file defined", 3000);
}

void STL2INP::on_pBtn_Output_2_pressed(){
	QString s = QFileDialog::getSaveFileName(
		this,
		"Choose a filename to save under",
		"E:/",
		"Data-File (*.dat)");
	if(s != ""){
		ui.pBtn_EditMesh->setEnabled(true);
		ui.lineEdit_Out_2->setText(s);
		ui.statusBar->showMessage("Output-File: Filename defined", 3000);
		on_pBtn_EditMesh_pressed();
	}
	else
		ui.statusBar->showMessage("Output-File: No file defined", 3000);
}

void STL2INP::on_lineEdit_Out_2_returnPressed(){
	QString s = ui.lineEdit_Out_2->text();
	if(s != ""){
		ui.pBtn_EditMesh->setEnabled(true);
		if(!s.endsWith(".dat",Qt::CaseSensitivity(false)))
		{
			s.append(".dat");
			ui.lineEdit_Out_2->setText(s);
		}
		if (QFile::exists(s))
			if(QMessageBox::question(this,
				tr("Overwrite File? -- Application Name"),
				tr("A file called %1 already exists."
				"Do you want to overwrite it?")
				.arg(s),
				tr("&Yes"), tr("&No"),
				QString(), 0, 1)){
					ui.lineEdit_Out_2->setText("");
					ui.pBtn_EditMesh->setEnabled(false);
			}
			else{
				ui.statusBar->showMessage("Output-File: overwrite existing File", 3000);
				ui.pBtn_EditMesh->setEnabled(true);
			}
		else{
			ui.statusBar->showMessage("Output-File: Filename defined", 3000);
			on_pBtn_EditMesh_pressed();
		}
	}
	else
		ui.statusBar->showMessage("Output-File: No file defined", 3000);
}

void STL2INP::on_pBtn_Convert_pressed(){
	if(ui.lineEdit_In->text() == "")
		QMessageBox::warning(this,"ERROR", "No Input-File defined!",
		QMessageBox::Ok, QMessageBox::NoButton, QMessageBox::NoButton);
	else if(ui.lineEdit_Out->text() == "" && ui.lineEdit_Out_2->text() == "")
		QMessageBox::warning(this,"ERROR", "No Output-File defined!",
		QMessageBox::Ok, QMessageBox::NoButton, QMessageBox::NoButton);  
	else
	{
		UbFileInputASCII *fileInput = new UbFileInputASCII(std::string(ui.lineEdit_In->text().toAscii()));
		GbTriangularMesh3D *mesh = GbTriangularMesh3DCreator::readMeshFromSTLFile(fileInput, "Cube");
		ui.statusBar->showMessage("Input-File was read", 3000);
		delete fileInput;
		cout<<mesh->toString()<<endl;
		if(ui.checkBox_writeAVS->isChecked()){
			if(ui.checkBox_Binary->isChecked()){
				UbFileOutputBinary *fileOutput_AVS = new UbFileOutputBinary(std::string(ui.lineEdit_Out->text().toAscii()));
				mesh->writeAVSMesh(fileOutput_AVS);
				delete fileOutput_AVS;
			}
			else{
				UbFileOutputASCII *fileOutput_AVS = new UbFileOutputASCII(std::string(ui.lineEdit_Out->text().toAscii()));
				mesh->writeAVSMesh(fileOutput_AVS);
				delete fileOutput_AVS;
			}
			ui.statusBar->showMessage("wrote AVS-Output-File");
		}
		if(ui.checkBox_writeUM->isChecked()){
			cout<<"MinX:"<<mesh->getX1Minimum()<<endl;
			cout<<"MaxX:"<<mesh->getX1Maximum()<<endl;
			cout<<"MinY:"<<mesh->getX2Minimum()<<endl;
			cout<<"MaxY:"<<mesh->getX2Maximum()<<endl;
			cout<<"MinZ:"<<mesh->getX3Minimum()<<endl;
			cout<<"MaxZ:"<<mesh->getX3Maximum()<<endl;
			ui.statusBar->showMessage("start Writing Uniform-Mesh-File");
			double minX = 0.0;
			double minY = 0.0;    
			double minZ = 0.0;

			CoordinateTransformation3D *trafo = new CoordinateTransformation3D(minX, minY, minZ, delta, delta, delta);

			UbTiming time;
			time.initTiming();
			time.startTiming();

			ui.statusBar->showMessage("start Building FacetteGrid", 3000);
			OctFacettenGrid2 *facettegrid = new OctFacettenGrid2("FacettenGrid", nx1, nx2, nx3, startLevel, stopLevel, mesh, trafo);
			ui.statusBar->showMessage("end Building FacetteGrid", 3000);

			UbFileOutputASCII out("E:/DATA/test.inp");
			facettegrid->writeCellsToAVS(&out);

			time.endTiming();
			cout<<"Dauer:"<<time.getDuration()<<endl;
			cout<<"Number of cells:"<<facettegrid->getNumberOfCells()<<endl;
			cout<<"after generation ..."<<endl<<endl;                                         
			double mydouble=0.0;

			time.initTiming(); 
			time.startTiming();
			ui.statusBar->showMessage("start writing", 3000);
			UbFileOutputASCII *fileOutput_UM = new UbFileOutputASCII(std::string(ui.lineEdit_Out_2->text().toAscii()));
			facettegrid->writeToUniformGridFile2(fileOutput_UM);
			delete fileOutput_UM;
			time.endTiming();
			cout<<"Dauer:"<<time.getDuration()<<endl;
			int number = (int)facettegrid->getCells()->size();
			delete trafo;
			delete mesh;
			delete facettegrid;
			cout<<"Ready!!!"<<endl;
			ui.statusBar->showMessage("wrote Unstructured-Mesh Output-File", 3000);
		}
	}
}

void STL2INP::on_checkBox_stateChanged(int)
{

}

void STL2INP::on_pBtn_EditMesh_pressed()
{
	QDefineUniformMesh *meshdef = new QDefineUniformMesh(this); 
	meshdef->setStartLevel(startLevel);
	meshdef->setStopLevel(stopLevel);
	meshdef->setDelta(delta);
	meshdef->setNX1(nx1);
	meshdef->setNX2(nx2);
	meshdef->setNX3(nx3);
	meshdef->exec();

	startLevel = meshdef->getStartLevel();
	stopLevel = meshdef->getStopLevel();
	//cout<<"Start-Level: "<<startLevel<<"  Stop-Level: "<<stopLevel<<endl;
	delta = meshdef->getDelta();
	//cout<<"Delta: "<<delta<<endl;
	nx1 = meshdef->getNX1();
	nx2 = meshdef->getNX2();
	nx3 = meshdef->getNX3();
	//cout<<"nx1: "<<nx1<<"  nx2: "<<nx2<<"  nx3: "<<nx3<<endl;
	delete meshdef;
}