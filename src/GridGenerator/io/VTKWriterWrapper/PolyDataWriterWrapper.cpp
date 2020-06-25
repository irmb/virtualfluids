//#include "PolyDataWriterWrapper.h"
//
//#include <VTKWriter/PolyDataWriter/PolyDataWriterImp.h>
//
//#include <GridGenerator/geometries/Arrow/Arrow.h>
//#include <GridGenerator/geometries/Vertex/Vertex.h>
//
//
//PolyDataWriterWrapper::PolyDataWriterWrapper()
//{
//	writer = std::shared_ptr<PolyDataWriter>(new PolyDataWriterImp());
//}
//
//PolyDataWriterWrapper::~PolyDataWriterWrapper()
//{
//
//}
//
//void PolyDataWriterWrapper::addVectorArrow(std::shared_ptr<const Arrow> arrow)
//{
//	double startPoint[] = { arrow->getStart()->x, arrow->getStart()->y, arrow->getStart()->z };
//	double endPoint[] = { arrow->getEnd()->x, arrow->getEnd()->y, arrow->getEnd()->z };
//	writer->addVectorArrow(startPoint, endPoint);
//}
//
//void PolyDataWriterWrapper::writePolyDataToFile(const std::string& filename) const
//{
//	writer->writePolyDataToFile(filename);
//}
