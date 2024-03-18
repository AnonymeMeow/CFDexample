cmake_minimum_required(VERSION 3.8)
project(SCARF)

aux_source_directory(${CMAKE_SOURCE_DIR} scarf_src)

link_libraries(m dl)

function(add_target target_name)
    add_executable(${target_name}_scarf ${scarf_src})

    target_compile_definitions(${target_name}_scarf PUBLIC ${target_name})
    set_target_properties(${target_name}_scarf PROPERTIES OUTPUT_NAME scarf)

    add_custom_target(${target_name}_build_and_run
        DEPENDS ${target_name}_scarf
        COMMAND ${CMAKE_SOURCE_DIR}/run.sh ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

add_subdirectory(Cavity)
add_subdirectory(Plate)
add_subdirectory(Tube)