
vf_add_library(BUILDTYPE shared DEPENDS basics GksMeshAdapter)

linkCUDA()
linkMPI()
linkOpenMP()
