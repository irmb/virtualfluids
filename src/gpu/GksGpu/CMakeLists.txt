
vf_add_library(BUILDTYPE shared PRIVATE_LINK basics GksMeshAdapter)

linkCUDA()
linkMPI()

vf_get_library_name(library_name)
linkOpenMP(${library_name})
