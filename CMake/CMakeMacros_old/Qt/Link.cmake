macro(linkQt targetName)

	find_package(Qt5Core)
	find_package(Qt5Widgets)
	find_package(Qt5Gui)
	find_package(Qt5PrintSupport)
	find_package(Qt5Charts)

	include_directories(${QT5Widgets_INCLUDES})
	include_directories(${QT5Core_INCLUDES})
	include_directories(${QT5Gui_INCLUDES})
	include_directories(${QT5PrintSupport_INCLUDES})

	add_definitions(${Qt5Widgets_DEFINITIONS})
	add_definitions(${Qt5Core_DEFINITIONS})
	add_definitions(${Qt5Gui_DEFINITIONS})
	add_definitions(${Qt5PrintSupport_DEFINITIONS})

	target_link_libraries(${targetName} Qt5::Widgets)
	target_link_libraries(${targetName} Qt5::Core)
	target_link_libraries(${targetName} Qt5::PrintSupport)
	target_link_libraries(${targetName} Qt5::Charts)

endmacro(linkQt)
