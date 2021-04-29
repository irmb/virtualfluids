//#include "Init/SetParameter.h"
//#include "Interface_OpenFOAM/Interface.h"
//
//////////////////////////////////////////////////////////////////////////////////
//void setParameters(Parameter* para, VF::GPU::Communicator* comm, std::string &cstr)
//{
//   ConfigFile cf(cstr.c_str());
//   if ( !cf.read() )
//   {
//      std::string exceptionText = "Unable to read configuration file\n";
//      throw exceptionText;
//   }
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   para->setMaxDev(              StringUtil::toInt(cf.getValue( "NumberOfDevices" )));
//   para->setMyID(                comm->getPID());                                     
//   para->setNumprocs(            comm->getNummberOfProcess());                        
//   para->setDevices(             StringUtil::toVector<int>(cf.getValue( "Devices" )));
//   devCheck(                     comm->mapCudaDevice(para->getMyID(), para->getNumprocs(), para->getDevices(), para->getMaxDev()));
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   std::string _path = cf.getValue( "Path" );
//   std::string _prefix = cf.getValue( "Prefix" );
//   std::string _gridpath = cf.getValue( "GridPath" );
//   std::string gridPath = getGridPath(para, _gridpath);
//   para->setOutputPath(          _path);
//   para->setOutputPrefix(        _prefix);
//   para->setFName(               _path + "/" + _prefix);
//   para->setPrintFiles(          false);                                                                  
//   para->setPrintFiles(          StringUtil::toBool( cf.getValue( "WriteGrid" )));                        
//   para->setGeometryValues(      StringUtil::toBool( cf.getValue( "GeometryValues" )));                        
//   para->setCalc2ndOrderMoments( StringUtil::toBool( cf.getValue( "calc2ndOrderMoments" )));                        
//   para->setCalc3rdOrderMoments( StringUtil::toBool( cf.getValue( "calc3rdOrderMoments" )));                        
//   para->setCalcHighOrderMoments(StringUtil::toBool( cf.getValue( "calcHigherOrderMoments" )));                        
//   para->setReadGeo(             StringUtil::toBool( cf.getValue( "ReadGeometry" )));                     
//   para->setCalcMedian(          StringUtil::toBool( cf.getValue( "calcMedian" )));                       
//   para->setConcFile(            StringUtil::toBool( cf.getValue( "UseConcFile" )));                       
//   para->setUseMeasurePoints(    StringUtil::toBool( cf.getValue( "UseMeasurePoints")));
//   para->setUseWale(             StringUtil::toBool( cf.getValue( "UseWale" )));
//   para->setSimulatePorousMedia( StringUtil::toBool( cf.getValue( "SimulatePorousMedia" )));
//   para->setD3Qxx(               StringUtil::toInt(  cf.getValue( "D3Qxx" )));
//   para->setMaxLevel(            StringUtil::toInt(  cf.getValue( "NOGL" )));                             
//   para->setTEnd(                StringUtil::toInt(  cf.getValue( "TimeEnd" )));                          
//   para->setTOut(                StringUtil::toInt(  cf.getValue( "TimeOut" )));                          
//   para->setTStartOut(			 StringUtil::toInt(  cf.getValue( "TimeStartOut")));
//   para->setTimeCalcMedStart(    StringUtil::toInt(  cf.getValue( "TimeStartCalcMedian" )));
//   para->setTimeCalcMedEnd(      StringUtil::toInt(  cf.getValue( "TimeEndCalcMedian" )));                  
//   para->setPressInID(           StringUtil::toInt(  cf.getValue( "PressInID" )));                          
//   para->setPressOutID(          StringUtil::toInt(  cf.getValue( "PressOutID" )));                         
//   para->setPressInZ(            StringUtil::toInt(  cf.getValue( "PressInZ" )));                           
//   para->setPressOutZ(           StringUtil::toInt(  cf.getValue( "PressOutZ" )));                          
//   //////////////////////////////////////////////////////////////////////////
//   para->setDiffOn(              StringUtil::toBool( cf.getValue( "DiffOn" )));                           
//   para->setDiffMod(             StringUtil::toInt(  cf.getValue( "DiffMod" )));                          
//   para->setDiffusivity(         StringUtil::toFloat(cf.getValue( "Diffusivity" )));                      
//   para->setTemperatureInit(     StringUtil::toFloat(cf.getValue( "Temp" )));                             
//   para->setTemperatureBC(       StringUtil::toFloat(cf.getValue( "TempBC" )));                           
//   //////////////////////////////////////////////////////////////////////////
//   para->setViscosity(           StringUtil::toFloat(cf.getValue( "Viscosity_LB" )));                     
//   para->setVelocity(            StringUtil::toFloat(cf.getValue( "Velocity_LB" )));                      
//   para->setViscosityRatio(      StringUtil::toFloat(cf.getValue( "Viscosity_Ratio_World_to_LB" )));      
//   para->setVelocityRatio(       StringUtil::toFloat(cf.getValue( "Velocity_Ratio_World_to_LB" )));       
//   para->setDensityRatio(        StringUtil::toFloat(cf.getValue( "Density_Ratio_World_to_LB" )));        
//   para->setPressRatio(          StringUtil::toFloat(cf.getValue( "Delta_Press" )));                      
//   para->setRealX(               StringUtil::toFloat(cf.getValue( "SliceRealX" )));                       
//   para->setRealY(               StringUtil::toFloat(cf.getValue( "SliceRealY" )));                       
//   para->setFactorPressBC(       StringUtil::toFloat(cf.getValue( "dfpbc" )));                      
//   para->setGeometryFileC(       cf.getValue( "GeometryC" ));                                             
//   para->setGeometryFileM(       cf.getValue( "GeometryM" ));                                             
//   para->setGeometryFileF(       cf.getValue( "GeometryF" ));                                             
//   //////////////////////////////////////////////////////////////////////////
//   para->setgeoVec(              gridPath + cf.getValue( "geoVec" ));
//   para->setcoordX(              gridPath + cf.getValue( "coordX" ));
//   para->setcoordY(              gridPath + cf.getValue( "coordY" ));
//   para->setcoordZ(              gridPath + cf.getValue( "coordZ" ));
//   para->setneighborX(           gridPath + cf.getValue( "neighborX" ));
//   para->setneighborY(           gridPath + cf.getValue( "neighborY" ));
//   para->setneighborZ(           gridPath + cf.getValue( "neighborZ" ));
//   para->setscaleCFC(            gridPath + cf.getValue( "scaleCFC" ));
//   para->setscaleCFF(            gridPath + cf.getValue( "scaleCFF" ));
//   para->setscaleFCC(            gridPath + cf.getValue( "scaleFCC" ));
//   para->setscaleFCF(            gridPath + cf.getValue( "scaleFCF" ));
//   para->setscaleOffsetCF(       gridPath + cf.getValue( "scaleOffsetCF" ));
//   para->setscaleOffsetFC(       gridPath + cf.getValue( "scaleOffsetFC" ));
//   para->setgeomBoundaryBcQs(    gridPath + cf.getValue( "geomBoundaryBcQs" ));
//   para->setgeomBoundaryBcValues(gridPath + cf.getValue( "geomBoundaryBcValues" ));
//   para->setinletBcQs(           gridPath + cf.getValue( "inletBcQs"        ));
//   para->setinletBcValues(       gridPath + cf.getValue( "inletBcValues"    ));
//   para->setoutletBcQs(          gridPath + cf.getValue( "outletBcQs"       ));
//   para->setoutletBcValues(      gridPath + cf.getValue( "outletBcValues"   ));
//   para->settopBcQs(             gridPath + cf.getValue( "topBcQs"          ));
//   para->settopBcValues(         gridPath + cf.getValue( "topBcValues"      ));
//   para->setbottomBcQs(          gridPath + cf.getValue( "bottomBcQs"       ));
//   para->setbottomBcValues(      gridPath + cf.getValue( "bottomBcValues"   ));
//   para->setfrontBcQs(           gridPath + cf.getValue( "frontBcQs"        ));
//   para->setfrontBcValues(       gridPath + cf.getValue( "frontBcValues"    ));
//   para->setbackBcQs(            gridPath + cf.getValue( "backBcQs"         ));
//   para->setbackBcValues(        gridPath + cf.getValue( "backBcValues"     ));
//   para->setnumberNodes(         gridPath + cf.getValue( "numberNodes"      ));
//   para->setLBMvsSI(             gridPath + cf.getValue( "LBMvsSI"          ));
//   //////////////////////////////gridPath + ////////////////////////////////////////////
//   para->setmeasurePoints(       gridPath + cf.getValue( "measurePoints" ));
//   para->setpropellerValues(	 gridPath + cf.getValue( "propellerValues"  ));
//   para->setclockCycleForMP(     StringUtil::toFloat(cf.getValue( "measureClockCycle" )));
//   para->settimestepForMP(       StringUtil::toInt(cf.getValue( "measureTimestep" )));
//   para->setcpTop(               gridPath + cf.getValue( "cpTop"            ));
//   para->setcpBottom(            gridPath + cf.getValue( "cpBottom"         ));
//   para->setcpBottom2(           gridPath + cf.getValue( "cpBottom2"        ));
//   para->setConcentration(       gridPath + cf.getValue( "Concentration"    ));
//   //////////////////////////////////////////////////////////////////////////
//   //Normals - Geometry
//   para->setgeomBoundaryNormalX(    gridPath + cf.getValue( "geomBoundaryNormalX" ));
//   para->setgeomBoundaryNormalY(    gridPath + cf.getValue( "geomBoundaryNormalY" ));
//   para->setgeomBoundaryNormalZ(    gridPath + cf.getValue( "geomBoundaryNormalZ" ));
//   //Normals - Inlet
//   para->setInflowBoundaryNormalX(    gridPath + cf.getValue( "inletBoundaryNormalX" ));
//   para->setInflowBoundaryNormalY(    gridPath + cf.getValue( "inletBoundaryNormalY" ));
//   para->setInflowBoundaryNormalZ(    gridPath + cf.getValue( "inletBoundaryNormalZ" ));
//   //Normals - Outlet
//   para->setOutflowBoundaryNormalX(    gridPath + cf.getValue( "outletBoundaryNormalX" ));
//   para->setOutflowBoundaryNormalY(    gridPath + cf.getValue( "outletBoundaryNormalY" ));
//   para->setOutflowBoundaryNormalZ(    gridPath + cf.getValue( "outletBoundaryNormalZ" ));
//   //////////////////////////////////////////////////////////////////////////
//   //Forcing
//   para->setForcing(StringUtil::toFloat(cf.getValue( "ForcingX")), StringUtil::toFloat(cf.getValue( "ForcingY")), StringUtil::toFloat(cf.getValue( "ForcingZ")));
//   //////////////////////////////////////////////////////////////////////////
//   //Particles
//   para->setCalcParticles(     StringUtil::toBool( cf.getValue( "calcParticles"     )));                             
//   para->setParticleBasicLevel(StringUtil::toInt(  cf.getValue( "baseLevel"         )));                             
//   para->setParticleInitLevel( StringUtil::toInt(  cf.getValue( "initLevel"         )));                             
//   para->setNumberOfParticles( StringUtil::toInt(  cf.getValue( "numberOfParticles" )));                             
//   para->setneighborWSB(                gridPath + cf.getValue( "neighborWSB"       ));
//   para->setStartXHotWall(     StringUtil::toDouble(cf.getValue("startXHotWall"     )));
//   para->setEndXHotWall(       StringUtil::toDouble(cf.getValue("endXHotWall"       )));
//   //////////////////////////////////////////////////////////////////////////
//   //for Multi GPU
//   if (para->getNumprocs()>1)
//   {
//	   ////////////////////////////////////////////////////////////////////////////
//	   ////1D domain decomposition
//	   //std::vector<std::string> sendProcNeighbors;
//	   //std::vector<std::string> recvProcNeighbors;
//	   //for (int i = 0; i<para->getNumprocs();i++)
//	   //{
//		  // sendProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "s.dat");
//		  // recvProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "r.dat");
//	   //}
//	   //para->setPossNeighborFiles(sendProcNeighbors, "send");
//	   //para->setPossNeighborFiles(recvProcNeighbors, "recv");
//	   //////////////////////////////////////////////////////////////////////////
//	   //3D domain decomposition
//	   std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
//	   std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
//	   for (int i = 0; i<para->getNumprocs();i++)
//	   {
//		   sendProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xs.dat");
//		   sendProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Ys.dat");
//		   sendProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zs.dat");
//		   recvProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xr.dat");
//		   recvProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Yr.dat");
//		   recvProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zr.dat");
//	   }
//	   para->setPossNeighborFilesX(sendProcNeighborsX, "send");
//	   para->setPossNeighborFilesY(sendProcNeighborsY, "send");
//	   para->setPossNeighborFilesZ(sendProcNeighborsZ, "send");
//	   para->setPossNeighborFilesX(recvProcNeighborsX, "recv");
//	   para->setPossNeighborFilesY(recvProcNeighborsY, "recv");
//	   para->setPossNeighborFilesZ(recvProcNeighborsZ, "recv");
//   }
//   //////////////////////////////////////////////////////////////////////////
//   //para->setkFull(             cf.getValue( "kFull" ));
//   //para->setgeoFull(           cf.getValue( "geoFull" ));
//   //para->setnoSlipBcPos(       cf.getValue( "noSlipBcPos" ));
//   //para->setnoSlipBcQs(          cf.getValue( "noSlipBcQs" ));
//   //para->setnoSlipBcValues(      cf.getValue( "noSlipBcValues" ));
//   //para->setnoSlipBcValue(     cf.getValue( "noSlipBcValue" ));
//   //para->setslipBcPos(         cf.getValue( "slipBcPos" ));
//   //para->setslipBcQs(          cf.getValue( "slipBcQs" ));
//   //para->setslipBcValue(       cf.getValue( "slipBcValue" ));
//   //para->setpressBcPos(        cf.getValue( "pressBcPos" ));
//   //para->setpressBcQs(           cf.getValue( "pressBcQs" ));
//   //para->setpressBcValues(       cf.getValue( "pressBcValues" ));
//   //para->setpressBcValue(      cf.getValue( "pressBcValue" ));
//   //para->setvelBcQs(             cf.getValue( "velBcQs" ));
//   //para->setvelBcValues(         cf.getValue( "velBcValues" ));
//   //para->setpropellerCylinder( cf.getValue( "propellerCylinder" ));
//   //para->setpropellerQs(		 cf.getValue( "propellerQs"      ));
//   //para->setwallBcQs(            cf.getValue( "wallBcQs"         ));
//   //para->setwallBcValues(        cf.getValue( "wallBcValues"     ));
//   //para->setperiodicBcQs(        cf.getValue( "periodicBcQs"     ));
//   //para->setperiodicBcValues(    cf.getValue( "periodicBcValues" ));
//   //cout << "Try this: " << para->getgeomBoundaryBcValues() << endl;
//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //Restart
//   para->setTimeDoCheckPoint(    StringUtil::toInt(  cf.getValue( "TimeDoCheckPoint" )));
//   para->setTimeDoRestart(       StringUtil::toInt(  cf.getValue( "TimeDoRestart" )));
//   para->setDoCheckPoint(        StringUtil::toBool( cf.getValue( "DoCheckPoint" )));
//   para->setDoRestart(           StringUtil::toBool( cf.getValue( "DoRestart" )));
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   para->setGridX(               StringUtil::toVector<int>(cf.getValue( "GridX" )));                      //GridX = StringUtil::toVector<int>(cf.getValue( "GridX" ));          
//   para->setGridY(               StringUtil::toVector<int>(cf.getValue( "GridY" )));                      //GridY = StringUtil::toVector<int>(cf.getValue( "GridY" ));          
//   para->setGridZ(               StringUtil::toVector<int>(cf.getValue( "GridZ" )));                      //GridZ = StringUtil::toVector<int>(cf.getValue( "GridZ" ));
//   para->setDistX(               StringUtil::toVector<int>(cf.getValue( "DistX" )));                      //DistX = StringUtil::toVector<int>(cf.getValue( "DistX" ));
//   para->setDistY(               StringUtil::toVector<int>(cf.getValue( "DistY" )));                      //DistY = StringUtil::toVector<int>(cf.getValue( "DistY" ));
//   para->setDistZ(               StringUtil::toVector<int>(cf.getValue( "DistZ" )));                      //DistZ = StringUtil::toVector<int>(cf.getValue( "DistZ" )); 
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   para->setNeedInterface(       StringUtil::toVector<bool>(cf.getValue( "NeedInterface" )));
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Interface *config = new Interface(false);
//   config->setDimensions(para);
//   config->setBoundingBox(para);
//   delete config;
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   para->initParameter();
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   para->setRe(para->getVelocity() * (real)1.0 / para->getViscosity());
//   para->setPhi((real) 0.0);
//   para->setlimitOfNodesForVTK(30000000); //max 30 Million nodes per VTK file
//   if (para->getDoRestart())
//   {
//	   para->setStartTurn(para->getTimeDoRestart());
//   } 
//   else
//   {
//	   para->setStartTurn((unsigned int)0); //100000
//   }
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////
//std::string getGridPath(Parameter* para, std::string Gridpath)
//{
//	if(para->getNumprocs()==1)
//	{
//		return Gridpath + "/";
//	}
//	else
//	{
//		return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
//	}
//
//}
//////////////////////////////////////////////////////////////////////////////////
//
