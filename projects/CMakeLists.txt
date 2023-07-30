cmake_minimum_required(VERSION 3.13)

project(tvb-cpp_projects)
set(CMAKE_PROJECT_NAME "tvb-cpp_projects")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_definitions(${tvb-cpp_DEFINITIONS})

add_executable(stroke_explore stroke/stroke_explore.cpp)
add_executable(anesthesia_explore sc1/anesthesia/anesthesia_explore.cpp sc1/anesthesia/zerlaut_gaba.h)
add_executable(propofol_explore_G sc1/anesthesia/propofol_explore_G.cpp sc1/anesthesia/zerlaut_gaba.h)
add_executable(zerlaut_figure paper_tvbcpp/zerlaut_figure.cpp)

find_package(MKL CONFIG)

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELDEBPROFILE  "-O2 -march=native -pg -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -g -march=native -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DEIGEN_NO_DEBUG -DNDEBUG")
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(stroke_explore Threads::Threads)
    target_link_libraries(anesthesia_explore Threads::Threads)
    target_link_libraries(propofol_explore_G Threads::Threads)
    target_link_libraries(zerlaut_figure Threads::Threads)

    if (MKL_FOUND)
        add_definitions( -DEIGEN_USE_MKL_ALL )
        include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(stroke_explore  $<LINK_ONLY:MKL::MKL>)
        target_link_libraries(anesthesia_explore  $<LINK_ONLY:MKL::MKL>)
        target_link_libraries(zerlaut_figure  $<LINK_ONLY:MKL::MKL>)
        target_link_libraries(propofol_explore_G  $<LINK_ONLY:MKL::MKL>)
    endif()

endif()

find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(stroke_explore Boost::program_options Boost::filesystem)
target_link_libraries(anesthesia_explore Boost::program_options Boost::filesystem)
target_link_libraries(propofol_explore_G Boost::program_options Boost::filesystem)

message(STATUS ${CMAKE_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/../src)

#find_package(Python3 3.8 EXACT COMPONENTS Development NumPy REQUIRED)
#include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
#target_link_libraries(stroke_explore Python3::Python Python3::NumPy)
#target_link_libraries(anesthesia_explore Python3::Python Python3::NumPy)
#target_link_libraries(propofol_explore_G Python3::Python Python3::NumPy)
#target_link_libraries(zerlaut_figure Python3::Python Python3::NumPy)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( stroke_explore  ${ZLIB_LIBRARIES} )
target_link_libraries( anesthesia_explore  ${ZLIB_LIBRARIES} )
target_link_libraries( propofol_explore_G  ${ZLIB_LIBRARIES} )
target_link_libraries( zerlaut_figure  ${ZLIB_LIBRARIES} )

get_filename_component(BUILD_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
set(TVB_ROOT_CPP_LIBRARYDIR ${CMAKE_SOURCE_DIR}/../src/${BUILD_NAME})
get_filename_component(TVB_ROOT_CPP_LIBRARYDIR ${TVB_ROOT_CPP_LIBRARYDIR} ABSOLUTE)
message(STATUS ${TVB_ROOT_CPP_LIBRARYDIR})

find_library(TVB_ROOT_CPP NAMES tvb-cpp PATHS "${TVB_ROOT_CPP_LIBRARYDIR}")
message(STATUS ${TVB_ROOT_CPP})

target_link_libraries( stroke_explore  ${TVB_ROOT_CPP} )
target_link_libraries( anesthesia_explore  ${TVB_ROOT_CPP} )
target_link_libraries( propofol_explore_G  ${TVB_ROOT_CPP} )
target_link_libraries( zerlaut_figure  ${TVB_ROOT_CPP} )
