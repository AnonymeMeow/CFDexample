cmake_minimum_required(VERSION 3.8)
project(CFDexample)

aux_source_directory(${CMAKE_SOURCE_DIR} scarf_src)

set(MPI_DIR "C:/Program Files (x86)/Microsoft SDKs/MPI")

function(add_target target_name)
    add_executable(${target_name}_grid ${CMAKE_SOURCE_DIR}/grid2D.c)

    target_compile_definitions(${target_name}_grid PUBLIC ${target_name})
    set_target_properties(${target_name}_grid PROPERTIES OUTPUT_NAME grid)

    add_executable(${target_name}_scarf ${scarf_src})

    target_compile_definitions(${target_name}_scarf
        PUBLIC ${target_name}
        PUBLIC SCARF
    )
    target_include_directories(${target_name}_scarf PUBLIC ${MPI_DIR}/Include)
    target_link_libraries(${target_name}_scarf PUBLIC ${MPI_DIR}/Lib/x64/msmpi.lib)
    set_target_properties(${target_name}_scarf PROPERTIES OUTPUT_NAME scarf)

    add_custom_target(${target_name}_build_and_run
        DEPENDS ${target_name}_grid ${target_name}_scarf
        COMMAND ${CMAKE_SOURCE_DIR}/run.bat ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

add_subdirectory(Cavity)
add_subdirectory(Plate)
add_subdirectory(Tube)