cmake_minimum_required(VERSION 3.29)
project(Project2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

if(DEFINED ENV{GUROBI_HOME})
    set(GUROBI_HOME "$ENV{GUROBI_HOME}")
    message(STATUS "GUROBI_HOME found: ${GUROBI_HOME}")
else()
    message(FATAL_ERROR "GUROBI_HOME not defined")
endif()

include_directories(${GUROBI_HOME}/include)

link_directories(${GUROBI_HOME}/lib)


if(WIN32)
    include_directories(${GUROBI_HOME}/include)
    link_directories(${GUROBI_HOME}/lib)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(GUROBI_CXX_LIBRARY gurobi_c++mdd2017)
    else()
        set(GUROBI_CXX_LIBRARY gurobi_c++md2017)
    endif()

    set(GUROBI_C_LIBRARY gurobi120)

else()
    # In Linux/macOS
    set(GUROBI_CXX_LIBRARY gurobi_c++)
    set(GUROBI_C_LIBRARY gurobi120)
endif()

add_executable(Project2
        main.cpp
        problem.cpp
        problem.h
        solution.cpp
        solution.h
        globals.cpp
        globals.h
        exact.cpp
        exact.h
        solution_exact.h
        solution_exact.cpp
)
target_include_directories(Project2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Project2 PRIVATE ${GUROBI_CXX_LIBRARY} ${GUROBI_C_LIBRARY})

if(UNIX AND NOT APPLE)
    set_target_properties(Project2 PROPERTIES
            INSTALL_RPATH "${GUROBI_HOME}/lib"
            BUILD_WITH_INSTALL_RPATH TRUE)
endif()