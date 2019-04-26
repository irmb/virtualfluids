#include "VTKAdapter.h"
#include "VTKInterface.h"

#include <vtkImageData.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

#include <vtkResampleWithDataSet.h>
#include <vtkPNGWriter.h>

#include <vtkGeometryFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/Definitions/MemoryAccessPattern.h"
#include "GksGpu/FlowStateData/FlowStateData.cuh"
#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"

vtkGridPtr getVtkUnstructuredOctGrid( SPtr<DataBase> dataBase, bool excludeGhostCells )
{
    vtkGridPtr grid = vtkGridPtr::New();
 
    vtkPointsPtr points = vtkPointsPtr::New();
 
    //////////////////////////////////////////////////////////////////////////

    for( uint nodeIdx = 0; nodeIdx < dataBase->numberOfNodes; nodeIdx++ ){
  
        Vec3& node = dataBase->nodeCoordinates[ nodeIdx ];

        points->InsertNextPoint( node.x, node.y, node.z );
    }

    grid->SetPoints( points );

    //////////////////////////////////////////////////////////////////////////

    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ ){

        if( dataBase->isGhostCell( cellIdx ) && excludeGhostCells ) continue;

        vtkIdListPtr idList = vtkIdListPtr::New();

        idList->SetNumberOfIds( 8 );

        idList->SetId( 0, dataBase->cellToNode[ cellIdx ][ 0 ] );
        idList->SetId( 1, dataBase->cellToNode[ cellIdx ][ 1 ] );
        idList->SetId( 2, dataBase->cellToNode[ cellIdx ][ 2 ] );
        idList->SetId( 3, dataBase->cellToNode[ cellIdx ][ 3 ] );
        idList->SetId( 4, dataBase->cellToNode[ cellIdx ][ 4 ] );
        idList->SetId( 5, dataBase->cellToNode[ cellIdx ][ 5 ] );
        idList->SetId( 6, dataBase->cellToNode[ cellIdx ][ 6 ] );
        idList->SetId( 7, dataBase->cellToNode[ cellIdx ][ 7 ] );

        grid->InsertNextCell( 12, idList );
    }

    //////////////////////////////////////////////////////////////////////////

    return grid;
}

void addScalarIntCellData( vtkGridPtr grid, 
                        uint numberOfCells, 
                        std::string name, 
                        std::function<int(uint)> getData )
{
    vtkIntArrayPtr data = vtkIntArrayPtr::New();

    data->SetNumberOfComponents( 1 );

    data->SetName( name.c_str() );

    for( uint cellIdx = 0; cellIdx < numberOfCells; cellIdx++ ){
        data->InsertNextValue( getData(cellIdx) );
    }

    grid->GetCellData()->AddArray( data );
}

void addScalarRealCellData( vtkGridPtr grid, 
                        uint numberOfCells, 
                        std::string name, 
                        std::function<real(uint)> getData )
{
    vtkDoubleArrayPtr data = vtkDoubleArrayPtr::New();

    data->SetNumberOfComponents( 1 );

    data->SetName( name.c_str() );

    for( uint cellIdx = 0; cellIdx < numberOfCells; cellIdx++ ){
        data->InsertNextValue( getData(cellIdx) );
    }

    grid->GetCellData()->AddArray( data );
}

void addVectorCellData( vtkGridPtr grid, 
                        uint numberOfCells, 
                        std::string name, 
                        std::function<Vec3(uint)> getData )
{
    vtkDoubleArrayPtr data = vtkDoubleArrayPtr::New();

    data->SetNumberOfComponents( 3 );

    data->SetName( name.c_str() );

    for( uint cellIdx = 0; cellIdx < numberOfCells; cellIdx++ ){
        Vec3 vec = getData(cellIdx);
        double tupel[3] = {vec.x, vec.y, vec.z};
        data->InsertNextTuple(tupel);
    }

    grid->GetCellData()->AddArray( data );
}

void addBaseData(vtkGridPtr grid, SPtr<DataBase> dataBase, Parameters parameters)
{
    addScalarIntCellData( grid, dataBase->numberOfCells, "CellIdx", [&] (uint cellIdx) {
        return cellIdx;
    } );

    addScalarRealCellData( grid, dataBase->numberOfCells, "rho", [&] (uint cellIdx) {
        return dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
    } );

    addScalarRealCellData( grid, dataBase->numberOfCells, "T", [&] (uint cellIdx) {
                    
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
        cons.rhoS_1 = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
        cons.rhoS_2 = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);
        
#ifdef USE_PASSIVE_SCALAR
        return getT(prim);
#else // USE_PASSIVE_SCALAR
        return 1.0 / prim.lambda;
#endif // USE_PASSIVE_SCALAR
    } );

    addScalarRealCellData( grid, dataBase->numberOfCells, "lambda", [&] (uint cellIdx) {
                    
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
        cons.rhoS_1 = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
        cons.rhoS_2 = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

        return prim.lambda;
    } );

    addScalarRealCellData( grid, dataBase->numberOfCells, "p", [&] (uint cellIdx) {
                    
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
        cons.rhoS_1 = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
        cons.rhoS_2 = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

        return 0.5 * prim.rho / prim.lambda;
    } );

    addScalarIntCellData( grid, dataBase->numberOfCells, "GhostCell", [&] (uint cellIdx) -> int {
        return dataBase->isGhostCell( cellIdx );
    } );

    addScalarIntCellData( grid, dataBase->numberOfCells, "Level", [&] (uint cellIdx) {
        return dataBase->getCellLevel(cellIdx);
    } );
            
    addVectorCellData( grid, dataBase->numberOfCells, "Velocity", [&] (uint cellIdx) {
                    
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];

        PrimitiveVariables prim = toPrimitiveVariables( cons, parameters.K );

        return Vec3( prim.U, prim.V, prim.W );
    } );

#ifdef USE_PASSIVE_SCALAR
	addScalarRealCellData( grid, dataBase->numberOfCells, "PassiveScalar_1", [&] (uint cellIdx) {
	    return dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ]
             / dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells)   ];
	} );

	addScalarRealCellData( grid, dataBase->numberOfCells, "PassiveScalar_2", [&] (uint cellIdx) {
	    return dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ]
             / dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells)   ];
	} );

	addScalarRealCellData( grid, dataBase->numberOfCells, "rhoE", [&] (uint cellIdx) {
	    return dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];
	} );
#endif // USE_PASSIVE_SCALAR

}

void writeVtkUnstructuredGrid( vtkGridPtr grid, int mode, std::string filename )
{
    vtkWriterPtr writer = vtkWriterPtr::New();

    writer->SetDataMode(mode);

    filename += ".";
    filename += writer->GetDefaultFileExtension();

    writer->SetFileName( filename.c_str() );

    writer->SetInputData( grid );

    writer->Write();
}

void VF_PUBLIC writeVtkParallelUnstructuredGridSummaryFile(vtkGridPtr grid, std::string filename, uint mpiWorldSize)
{
    uint numberOfArrays = grid->GetCellData()->GetNumberOfArrays();

    const auto filenameWithoutPath=filename.substr( filename.find_last_of('/') + 1 );

    std::ofstream file;

    file.open( filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PCellData>" << std::endl;

    for( uint i = 0; i < numberOfArrays; i++ )
    {
        int typeID( grid->GetCellData()->GetArray(i)->GetDataType() );
        std::string name( grid->GetCellData()->GetArray(i)->GetName() );

        uint numberOfComponents = grid->GetCellData()->GetArray(i)->GetNumberOfComponents();

        std::string type;
        if( typeID == VTK_INT    ) type = "Int32";
        if( typeID == VTK_FLOAT  ) type = "Float32";
        if( typeID == VTK_DOUBLE ) type = "Float64";

        file << "      <PDataArray type=\"" << type << "\" Name=\"" << name << "\" NumberOfComponents=\"" << numberOfComponents << "\"/>" << std::endl;
    }

    file << "    </PCellData>" << std::endl;

    file << "    <PPoints>" << std::endl;
    file << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << std::endl;
    file << "    </PPoints>" << std::endl;

    for( uint rank = 0; rank < mpiWorldSize; rank++ )
    {
        file << "    <Piece Source=\"" << filenameWithoutPath << "_rank_" << rank << ".vtu\"/>" << std::endl;
    }

    file << "  </PUnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;

    //////////////////////////////////////////////////////////////////////////
}

rgbColor colorMapCoolToWarmExtended( double value, double min, double max )
{    
    // Color map exported from Paraview
    const double colorMap[36][3] = 
    /*  0 */  { { 0,                      0,                      0.34902              },
    /*  1 */    { 0.039216000000000001,   0.062744999999999995,   0.38039200000000001  },
    /*  2 */    { 0.062744999999999995,   0.117647,               0.41176499999999999  },
    /*  3 */    { 0.090195999999999998,   0.18431400000000001,    0.45097999999999999  },
    /*  4 */    { 0.12548999999999999,    0.26274500000000001,    0.50196099999999999  },
    /*  5 */    { 0.16078400000000001,    0.33725500000000003,    0.54117599999999999  },
    /*  6 */    { 0.20000000000000001,    0.39607799999999999,    0.56862699999999999  },
    /*  7 */    { 0.23921600000000001,    0.45490199999999997,    0.59999999999999998  },
    /*  8 */    { 0.286275,               0.52156899999999995,    0.65098              },
    /*  9 */    { 0.33725500000000003,    0.59215700000000004,    0.70196099999999995  },
    /* 10 */    { 0.388235,               0.65490199999999998,    0.74902000000000002  },
    /* 11 */    { 0.466667,               0.73725499999999999,    0.819608             },
    /* 12 */    { 0.57254899999999997,    0.819608,               0.87843099999999996  },
    /* 13 */    { 0.65490199999999998,    0.86666699999999997,    0.90980399999999995  },
    /* 14 */    { 0.75294099999999997,    0.91764699999999999,    0.94117600000000001  },
    /* 15 */    { 0.82352899999999996,    0.95686300000000002,    0.96862700000000002  },
    ///* 15 */    { 1.0,                    1.0,                    1.0                  },
    /* 16 */    { 0.98823499999999997,    0.96078399999999997,    0.90196100000000001  },
    ///* 16 */    { 1.0,                    1.0,                    1.0                  },

    ///* 17 */    { 1.0,                    1.0,                    1.0                  },
    /* 17 */    { 0.94117600000000001,    0.98431400000000002,    0.98823499999999997  },
    ///* 18 */    { 1.0,                    1.0,                    1.0                  },
    /* 18 */    { 0.98823499999999997,    0.94509799999999999,    0.85097999999999996  },
    ///* 19 */    { 1.0,                    1.0,                    1.0                  },
    /* 19 */    { 0.98039200000000004,    0.89803900000000003,    0.78431399999999996  },
    /* 20 */    { 0.96862700000000002,    0.83529399999999998,    0.69803899999999997  },
    /* 21 */    { 0.94901999999999997,    0.73333300000000001,    0.58823499999999995  },
    /* 22 */    { 0.92941200000000002,    0.65098,                0.50980400000000003  },
    /* 23 */    { 0.90980399999999995,    0.56470600000000004,    0.43529400000000001  },
    /* 24 */    { 0.87843099999999996,    0.45882400000000001,    0.352941             },
    /* 25 */    { 0.83921599999999996,    0.388235,               0.286275             },
    /* 26 */    { 0.76078400000000002,    0.29411799999999999,    0.21176500000000001  },
    /* 27 */    { 0.70196099999999995,    0.21176500000000001,    0.168627             },
    /* 28 */    { 0.65098,                0.156863,               0.129412             },
    /* 29 */    { 0.59999999999999998,    0.094117999999999993,   0.094117999999999993 },
    /* 30 */    { 0.54901999999999995,    0.066667000000000004,   0.098039000000000001 },
    /* 31 */    { 0.50196099999999999,    0.050979999999999998,   0.12548999999999999  },
    /* 32 */    { 0.45097999999999999,    0.054901999999999999,   0.17254900000000001  },
    /* 33 */    { 0.40000000000000002,    0.054901999999999999,   0.19215699999999999  },
    /* 34 */    { 0.34902,                0.070587999999999998,   0.21176500000000001  },
                { 0.34902,                0.070587999999999998,   0.21176500000000001  } };
    
    if( value < min )
        value = 0.0;
    else if ( value > max )
        value = 1.0;
    else
        value = ( value - min ) / ( max - min );

    unsigned int idx           = value * 34;
    double       interpolation = value * 34.0 - double( idx );

    rgbColor color;

    color.r = ( ( 1.0 - interpolation ) * colorMap[idx  ][0]
              +         interpolation   * colorMap[idx+1][0] ) * 200.0 * 1.2;

    color.g = ( ( 1.0 - interpolation ) * colorMap[idx  ][1]
              +         interpolation   * colorMap[idx+1][1] ) * 197.0 * 1.2;

    color.b = ( ( 1.0 - interpolation ) * colorMap[idx  ][2]
              +         interpolation   * colorMap[idx+1][2] ) * 189.0 * 1.2;

    return color;
}

void writePNG( vtkDataObject * inputData, int nx, int ny, double L, double H, std::string filename )
{
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions( nx, ny, 1 );
    image->SetSpacing( L / ( nx - 1 ), H / ( ny - 1 ), 0 );

    vtkSmartPointer<vtkResampleWithDataSet> resample = vtkSmartPointer<vtkResampleWithDataSet>::New();
    resample->SetSourceData( inputData );
    resample->SetInputData( image );
    resample->Update();

    vtkSmartPointer<vtkImageData> image2 = (vtkImageData*) resample->GetOutput();

    image2->GetPointData()->SetScalars( image2->GetPointData()->GetArray( 0 ) );

    vtkSmartPointer<vtkPNGWriter> writerPNG = vtkSmartPointer<vtkPNGWriter>::New();
    writerPNG->SetFileName( ( filename + ".png" ).c_str() );
    writerPNG->SetInputData( image2 );
    writerPNG->Write();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void writeVtkXML(std::shared_ptr<DataBase> dataBase, 
                 Parameters parameters, 
                 int mode, 
                 std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << filename << ".vtu" << " ... \n";

    vtkGridPtr grid = getVtkUnstructuredOctGrid(dataBase);

    addBaseData( grid, dataBase, parameters );

    writeVtkUnstructuredGrid( grid, vtkXMLWriter::Binary, filename );

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void VF_PUBLIC writeVtkXMLParallelSummaryFile(std::shared_ptr<DataBase> dataBase, Parameters parameters, std::string filename, uint mpiWorldSize)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << filename << ".pvtu" << " ... \n";

    vtkGridPtr grid = getVtkUnstructuredOctGrid(dataBase);

    addBaseData( grid, dataBase, parameters );

    writeVtkParallelUnstructuredGridSummaryFile( grid, filename, mpiWorldSize );

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void writeTurbulenceVtkXML(std::shared_ptr<DataBase> dataBase, 
                           std::shared_ptr<TurbulenceAnalyzer> turbulenceAnalyzer,
                           int mode, 
                           std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << filename << ".vtu" << " ... \n";

    vtkGridPtr grid = getVtkUnstructuredOctGrid(dataBase);

    addScalarIntCellData( grid, dataBase->numberOfCells, "CellIdx", [&] (uint cellIdx) {
        return cellIdx;
    } );

    addScalarIntCellData( grid, dataBase->numberOfCells, "GhostCell", [&] (uint cellIdx) -> int {
        return dataBase->isGhostCell( cellIdx );
    } );

    //////////////////////////////////////////////////////////////////////////

    addScalarRealCellData(grid, dataBase->numberOfCells, "U", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_U[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "V", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_V[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "W", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_W[ cellIdx ];
    });

    //////////////////////////////////////////////////////////////////////////

    addScalarRealCellData(grid, dataBase->numberOfCells, "UU", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_UU[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "VV", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_VV[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "WW", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_WW[ cellIdx ];
    });

    //////////////////////////////////////////////////////////////////////////

    addScalarRealCellData(grid, dataBase->numberOfCells, "UV", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_UV[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "UW", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_UW[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "VW", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_VW[ cellIdx ];
    });

    //////////////////////////////////////////////////////////////////////////

    addScalarRealCellData(grid, dataBase->numberOfCells, "T", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_T[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "TT", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_TT[ cellIdx ];
    });

    addScalarRealCellData(grid, dataBase->numberOfCells, "p", [&](uint cellIdx) {
        return turbulenceAnalyzer->h_p[ cellIdx ];
    });

    //////////////////////////////////////////////////////////////////////////

    writeVtkUnstructuredGrid( grid, vtkXMLWriter::Binary, filename );

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void VF_PUBLIC writeTurbulenceVtkXMLParallelSummaryFile(std::shared_ptr<DataBase> dataBase, std::shared_ptr<TurbulenceAnalyzer> turbulenceAnalyzer,Parameters parameters, std::string filename, uint mpiWorldSize)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << filename << ".pvtu" << " ... \n";

    vtkGridPtr grid = getVtkUnstructuredOctGrid(dataBase);
    
    //////////////////////////////////////////////////////////////////////////

    const auto filenameWithoutPath=filename.substr( filename.find_last_of('/') + 1 );

    std::ofstream file;

    file.open( filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PCellData>" << std::endl;

        file << "      <PDataArray type=\"" << "Int32"   << "\" Name=\"" << "CellIdx"   << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Int32"   << "\" Name=\"" << "GhostCell" << "\" NumberOfComponents=\"1\"/>" << std::endl;

        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "U"         << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "V"         << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "W"         << "\" NumberOfComponents=\"1\"/>" << std::endl;

        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "UU"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "VV"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "WW"        << "\" NumberOfComponents=\"1\"/>" << std::endl;

        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "UV"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "UW"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "VW"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "T"         << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "TT"        << "\" NumberOfComponents=\"1\"/>" << std::endl;
        file << "      <PDataArray type=\"" << "Float64" << "\" Name=\"" << "p"         << "\" NumberOfComponents=\"1\"/>" << std::endl;

    file << "    </PCellData>" << std::endl;

    file << "    <PPoints>" << std::endl;
    file << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << std::endl;
    file << "    </PPoints>" << std::endl;

    for( uint rank = 0; rank < mpiWorldSize; rank++ )
    {
        file << "    <Piece Source=\"" << filenameWithoutPath << "_rank_" << rank << ".vtu\"/>" << std::endl;
    }

    file << "  </PUnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void mapFlowField(std::shared_ptr<DataBase> base, std::shared_ptr<DataBase> target)
{
    vtkGridPtr gridBase   = getVtkUnstructuredOctGrid(base,   true);
    vtkGridPtr gridTarget = getVtkUnstructuredOctGrid(target, true);

    //////////////////////////////////////////////////////////////////////////

    vtkSmartPointer<vtkDoubleArray> rho  = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> rhoU = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> rhoE = vtkSmartPointer<vtkDoubleArray>::New();
            
    rho->SetNumberOfComponents ( 1 );
    rhoU->SetNumberOfComponents( 3 );
    rhoE->SetNumberOfComponents( 1 );

    rho->SetName ( "rho"  );
    rhoU->SetName( "rhoU" );
    rhoE->SetName( "rhoW" );

    for( uint cellIdx = 0; cellIdx < base->numberOfCells; cellIdx++ ){

        if( base->isGhostCell( cellIdx ) ) continue;
                    
        ConservedVariables cons;

        cons.rho  = base->dataHost[ RHO__(cellIdx, base->numberOfCells) ];
        cons.rhoU = base->dataHost[ RHO_U(cellIdx, base->numberOfCells) ];
        cons.rhoV = base->dataHost[ RHO_V(cellIdx, base->numberOfCells) ];
        cons.rhoW = base->dataHost[ RHO_W(cellIdx, base->numberOfCells) ];
        cons.rhoE = base->dataHost[ RHO_E(cellIdx, base->numberOfCells) ];

        rho->InsertNextTuple1 ( cons.rho );
        rhoU->InsertNextTuple3( cons.rhoU, cons.rhoV, cons.rhoW );
        rhoE->InsertNextTuple1( cons.rhoE );
    }

    gridBase->GetCellData()->AddArray( rho  );
    gridBase->GetCellData()->AddArray( rhoU );
    gridBase->GetCellData()->AddArray( rhoE );
        
#ifdef USE_PASSIVE_SCALAR

        vtkSmartPointer<vtkDoubleArray> dataS_1 = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> dataS_2 = vtkSmartPointer<vtkDoubleArray>::New();

        dataS_1->SetNumberOfComponents(1);
        dataS_2->SetNumberOfComponents(1);

        dataS_1->SetName("rhoS_1");
        dataS_2->SetName("rhoS_2");

        for (uint cellIdx = 0; cellIdx < base->numberOfCells; cellIdx++) {

            if (base->isGhostCell(cellIdx)) continue;

            dataS_1->InsertNextTuple1(base->dataHost[RHO_S_1(cellIdx, base->numberOfCells)]);
            dataS_2->InsertNextTuple1(base->dataHost[RHO_S_2(cellIdx, base->numberOfCells)]);
        }

        gridBase->GetCellData()->AddArray(dataS_1);
        gridBase->GetCellData()->AddArray(dataS_2);

#endif // USE_PASSIVE_SCALAR

    //////////////////////////////////////////////////////////////////////////

    vtkSmartPointer<vtkCellDataToPointData> cellDataToPointDataBase = vtkSmartPointer<vtkCellDataToPointData>::New();
    cellDataToPointDataBase->SetInputData( gridBase );
    cellDataToPointDataBase->Update();

    vtkSmartPointer<vtkGeometryFilter> gridToPolyDataBase = vtkSmartPointer<vtkGeometryFilter>::New();
    gridToPolyDataBase->SetInputConnection( cellDataToPointDataBase->GetOutputPort() );
    gridToPolyDataBase->Update();

    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputConnection( gridToPolyDataBase->GetOutputPort() );
    cleanPolyData->Update();

    vtkSmartPointer<vtkResampleWithDataSet> resampleWithDataSet = vtkSmartPointer<vtkResampleWithDataSet>::New();
    resampleWithDataSet->SetSourceConnection( cleanPolyData->GetOutputPort() );
    resampleWithDataSet->SetInputData(  gridTarget );
    resampleWithDataSet->Update();

    vtkSmartPointer<vtkPointDataToCellData> pointDataToCellDataTarget = vtkSmartPointer<vtkPointDataToCellData>::New();
    pointDataToCellDataTarget->SetInputConnection( resampleWithDataSet->GetOutputPort() );
    pointDataToCellDataTarget->Update();

    gridTarget = (vtkUnstructuredGrid*) pointDataToCellDataTarget->GetOutput();

    //////////////////////////////////////////////////////////////////////////

    for( uint cellIdx = 0, gridCellIdx = 0; cellIdx < target->numberOfCells; cellIdx++ ){

        if( target->isGhostCell( cellIdx ) ) continue;

        double  rho  = gridTarget->GetCellData()->GetArray(1)->GetTuple1(gridCellIdx);
        double* rhoU = gridTarget->GetCellData()->GetArray(2)->GetTuple3(gridCellIdx);
        double  rhoE = gridTarget->GetCellData()->GetArray(3)->GetTuple1(gridCellIdx);

        target->dataHost[ RHO__(cellIdx, target->numberOfCells) ] = rho;
        target->dataHost[ RHO_U(cellIdx, target->numberOfCells) ] = rhoU[0];
        target->dataHost[ RHO_V(cellIdx, target->numberOfCells) ] = rhoU[1];
        target->dataHost[ RHO_W(cellIdx, target->numberOfCells) ] = rhoU[2];
        target->dataHost[ RHO_E(cellIdx, target->numberOfCells) ] = rhoE;

#ifdef USE_PASSIVE_SCALAR
        {
            double  rhoS_1 = gridTarget->GetCellData()->GetArray(4)->GetTuple1(gridCellIdx);
            double  rhoS_2 = gridTarget->GetCellData()->GetArray(5)->GetTuple1(gridCellIdx);

            target->dataHost[RHO_S_1(cellIdx, target->numberOfCells)] = rhoS_1;
            target->dataHost[RHO_S_2(cellIdx, target->numberOfCells)] = rhoS_2;
        }
#endif // USE_PASSIVE_SCALAR

        gridCellIdx++;
    }
}