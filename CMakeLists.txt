cmake_minimum_required(VERSION 3.8)
project(CFDexample)

add_executable(grid "grid2D.c")

add_executable(scarf
    "bound.c"
    "chemdata.h"
    "comm.h"
    "end.c"
    "flux.c"
    "jobbody.c"
    "main.c"
    "read.c"
    "save.c"
    "set.c"
    "start.c"
    "thermochem.c"
    "tmarch.c"
    "vflux.c"
)

set(MPI_DIR "C:/Program Files (x86)/Microsoft SDKs/MPI")

target_include_directories(scarf PUBLIC ${MPI_DIR}/Include)

target_link_libraries(scarf PUBLIC ${MPI_DIR}/Lib/x64/msmpi.lib)

function(add_target target_name)
    target_compile_definitions(scarf PUBLIC ${target_name})
    add_custom_target(${target_name}_build_and_run
        DEPENDS grid scarf
        COMMAND ${CMAKE_SOURCE_DIR}/run.bat ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}
    )
endfunction()

add_subdirectory(Plate)