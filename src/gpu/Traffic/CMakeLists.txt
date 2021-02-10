project(Traffic LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE shared PRIVATE_LINK GridGenerator basics)
