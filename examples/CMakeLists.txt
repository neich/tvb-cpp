cmake_minimum_required(VERSION 3.13)

project(tvb_cpp_examples)
set(CMAKE_PROJECT_NAME "tvb-root-cpp-examples")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(test_simulation test_simulation.cpp)
add_executable(parameter_sweep parameter_sweep.cpp)

find_package(MKL CONFIG)

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELDEBPROFILE  "-O2 -march=native -pg -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DEIGEN_NO_DEBUG -DNDEBUG")
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(test_simulation Threads::Threads)
    target_link_libraries(parameter_sweep Threads::Threads)

    if (MKL_FOUND)
        add_definitions( -DEIGEN_USE_MKL_ALL )
        include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(test_simulation  $<LINK_ONLY:MKL::MKL>)
        target_link_libraries(parameter_sweep  $<LINK_ONLY:MKL::MKL>)
    endif()

endif()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:preprocessor")
endif()


find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(test_simulation Boost::program_options Boost::filesystem)
target_link_libraries(parameter_sweep Boost::program_options  Boost::filesystem)

include_directories(${CMAKE_SOURCE_DIR}/../src)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( test_simulation  ${ZLIB_LIBRARIES} )
target_link_libraries( parameter_sweep ${ZLIB_LIBRARIES} )

get_filename_component(BUILD_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
set(TVB_ROOT_CPP_LIBRARYDIR ${CMAKE_SOURCE_DIR}/../src/${BUILD_NAME})
get_filename_component(TVB_ROOT_CPP_LIBRARYDIR ${TVB_ROOT_CPP_LIBRARYDIR} ABSOLUTE)
find_library(TVB_ROOT_CPP NAMES tvb-cpp PATHS "${TVB_ROOT_CPP_LIBRARYDIR}")

target_link_libraries( test_simulation  ${TVB_ROOT_CPP} )
target_link_libraries( parameter_sweep ${TVB_ROOT_CPP} )
