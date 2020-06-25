macro(bundleTargetForMacOSX targetName)

    set_target_properties(${targetName} PROPERTIES MACOSX_BUNDLE ON)
    set_target_properties(${targetName} PROPERTIES MACOSX_BUNDLE_INFO_PLIST  ${CMAKE_SOURCE_DIR}/cmake/MacOSXBundle/Info.plist)

endmacro(bundleTargetForMacOSX)