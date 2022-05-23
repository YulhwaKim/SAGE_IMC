cmake_minimum_required(VERSION 2.8)
project(cim_dse)

find_package(OpenMP REQUIRED)

option(CMAKE_EXPORT_COMPILE_COMMANDS "" 1)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("--std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -std=c++11") # -w disables warnings
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++ support. Please use a different C++ compiler.")
endif()

include_directories("${CMAKE_SOURCE_DIR}/CircuitModule/include")
include_directories("${CMAKE_SOURCE_DIR}/ArchitectureTemplate/include")
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/Scheduler/include")

add_subdirectory("src")
add_subdirectory("CircuitModule/src")
add_subdirectory("ArchitectureTemplate/src")
add_subdirectory("Scheduler/src")

# set custom library
set(CM_LIBS circuit_module)
set(HT_LIBS architecture_template)
set(SC_LIBS scheduler)
set(TOP_LIBS top_source ${OpenMP_LIBRARIES})

# merge libraray
set(DSE_LIBS ${CMAKE_BINARY_DIR}/libdse/libdse.a)
add_custom_target(combined
        COMMAND ar -x $<TARGET_FILE:circuit_module>
        COMMAND ar -x $<TARGET_FILE:architecture_template>
        COMMAND ar -x $<TARGET_FILE:scheduler>
        COMMAND ar -x $<TARGET_FILE:top_source>
        COMMAND ar -qcs ${DSE_LIBS} *.o
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/libdse
        DEPENDS circuit_module architecture_template scheduler top_source
    )

add_library(dse STATIC IMPORTED GLOBAL)
add_dependencies(dse combined)

set_target_properties(dse
        PROPERTIES
        IMPORTED_LOCATION ${DSE_LIBS}
    )


# make runable file
add_executable(arch_generator_for_layer arch_generator_for_layer.cpp)
target_link_libraries(arch_generator_for_layer dse ${OpenMP_LIBRARIES})

add_executable(arch_generator_for_net arch_generator_for_net.cpp)
target_link_libraries(arch_generator_for_net dse ${OpenMP_LIBRARIES})

add_executable(main main.cpp)
target_link_libraries(main dse ${OpenMP_LIBRARIES})

add_executable(main_iter main_iter.cpp)
target_link_libraries(main_iter dse ${OpenMP_LIBRARIES})

