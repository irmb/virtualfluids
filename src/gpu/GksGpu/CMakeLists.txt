
vf_add_library(BUILDTYPE shared PRIVATE_LINK basics GksMeshAdapter)

linkCUDA()
linkMPI()
linkOpenMP()
