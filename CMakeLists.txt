cmake_minimum_required(VERSION 3.8)
project(CFDexample)

add_executable(grid "grid2D.c")

aux_source_directory(${CMAKE_SOURCE_DIR} scarf_src)

set(MPI_DIR "C:/Program Files (x86)/Microsoft SDKs/MPI")

include_directories(${MPI_DIR}/Include)

link_libraries(${MPI_DIR}/Lib/x64/msmpi.lib)

function(add_target target_name)
    add_executable(${target_name}_scarf ${scarf_src})
    target_compile_definitions(${target_name}_scarf PUBLIC ${target_name})
    set_target_properties(${target_name}_scarf PROPERTIES OUTPUT_NAME scarf)
    add_custom_target(${target_name}_build_and_run
        DEPENDS grid ${target_name}_scarf
        COMMAND ${CMAKE_SOURCE_DIR}/run.bat ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

add_subdirectory(Cavity)
add_subdirectory(Plate)
add_subdirectory(Tube)