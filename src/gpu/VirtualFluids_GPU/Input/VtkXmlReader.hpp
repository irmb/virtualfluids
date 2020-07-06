/* 
 * File:   VtkXmlReader.hpp
 * Author: kucher
 *
 * Created on 3. August 2010, 15:59
 *
 * extract PointData from .VTK file
 */

#ifndef VTK_XML_READER_H
#define VTK_XML_READER_H

#include "Utilities/StringUtil.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iostream>

namespace vtk_xml_reader
{
	template <class T>
      
  // T* getPointData(const std::string &filename, const std::string &dataArrayName, const std::string &gridTyp){
  /*std::vector<T>*/ 
   void getPointData(const std::string &filename, const std::string &dataArrayName, const std::string &gridTyp, T* pointData, const int &pointDataSize){

		//std::vector<T> pointData;
		std::vector<std::string> pointDataStrings;
		std::string str;

		// Create empty property tree object
		using boost::property_tree::ptree;
		ptree pt;


		try
		{
			// Load XML file and put its contents in property tree.
			// No namespace qualification is needed, because of Koenig
			// lookup on the second argument. If reading fails, exception
			// is thrown.
			read_xml(filename, pt);

         std::string required_tag = "VTKFile." + gridTyp + ".Piece.PointData";

			//BOOST_FOREACH(ptree::value_type &v, pt.get_child("VTKFile.UnstructuredGrid.Piece.PointData")){
         BOOST_FOREACH(ptree::value_type &v, pt.get_child(required_tag)){

				if(v.first == "DataArray"){
					if(v.second.get<std::string>("<xmlattr>.Name") == dataArrayName)
						str = v.second.data();
				}
			}
		}
		catch (boost::property_tree::ptree_error  e)
		{
			std::cout << e.what() << std::endl;
			exit(0);
		}

      boost::algorithm::split(pointDataStrings, str, boost::is_any_of("\t\n "));

      int i = 0;
		BOOST_FOREACH(std::string s, pointDataStrings){
			if (s != "")
         {
				//pointData.push_back(StringUtil::fromString<T>(s));
            pointData[i] = StringUtil::fromString<T>(s);
            i++;
            if (i > pointDataSize)
            {
               std::string exceptionText = "Bad allocation: pointDataStrings > pointDataSize";
               throw exceptionText;
            }
         }
		}

		//return &pointData[0];
      //return pointData;
      
	}
}

#endif	

