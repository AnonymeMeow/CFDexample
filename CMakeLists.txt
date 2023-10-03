cmake_minimum_required (VERSION 3.8)
project (CFD)

set (CMAKE_CXX_STANDARD 17)

aux_source_directory (./ SOURCE_FILES)
add_executable (Cavity ${SOURCE_FILES})