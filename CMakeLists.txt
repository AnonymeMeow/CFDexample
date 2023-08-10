cmake_minimum_required (VERSION 3.8)
project (CFD)

set (CMAKE_CXX_STANDARD 17)

include_directories ("C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
link_directories ("C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64")
aux_source_directory (./ SOURCE_FILES)
add_executable (Cavity ${SOURCE_FILES})
target_link_libraries (Cavity msmpi)