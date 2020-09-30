project(NumericalTestPostProcessing)

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluids_GPU basics cudart_static)

vf_get_library_name(library_name)
groupTarget(${library_name} ${testFolder})
