project(NumericalTestPostProcessing)

vf_add_library(BUILDTYPE binary PRIVATE_LINK gpu_core basics cudart_static)

vf_get_library_name(library_name)
groupTarget(${library_name} ${testFolder})
