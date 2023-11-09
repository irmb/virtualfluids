PROJECT(MusselOyster LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES MusselOyster.cpp)
