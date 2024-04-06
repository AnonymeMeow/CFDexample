cmake_minimum_required(VERSION 3.8)
project(SCARF)

aux_source_directory(${CMAKE_SOURCE_DIR} scarf_src)

link_libraries(m dl)

add_executable(scarf ${scarf_src})

add_custom_target(build_and_run
    DEPENDS scarf
    COMMAND ${CMAKE_SOURCE_DIR}/run.sh ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)