cmake_minimum_required(VERSION 3.13)

project(tvb_cpp_examples)
set(CMAKE_PROJECT_NAME "tvb-root-cpp-examples")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (MSVC)
#    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "/Od" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELINTODEB  "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FS")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} /Ox")
endif()

set(CMAKE_CXX_FLAGS_RELINTODEB  "")

add_executable(fic_explore fic_explore.cpp)
add_executable(parameter_sweep parameter_sweep.cpp)
add_executable(test_simulation test_simulation.cpp)
add_executable(stroke_explore stroke/stroke_explore.cpp)

find_package(MKL CONFIG)

if (MKL_FOUND)
    add_definitions( -DEIGEN_USE_MKL_ALL )
    include_directories(${MKL_H})
endif()

find_package(BLAS)

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g -DEIGEN_NO_DEBUG" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -ffast-math -DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
    # set(CMAKE_CXX_FLAGS_DEBUG  "-g -save-temps" CACHE STRING "" FORCE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(fic_explore Threads::Threads)
    target_link_libraries(stroke_explore Threads::Threads)
    target_link_libraries(parameter_sweep Threads::Threads)

    if (MKL_FOUND)
        target_link_libraries(fic_explore mkl_rt)
        target_link_libraries(parameter_sweep mkl_rt)
        target_link_libraries(test_simulation mkl_rt)
        target_link_libraries(stroke_explore mkl_rt)
    endif()
    if (BLAS_FOUND)
        target_link_libraries(fic_explore mkl_rt)
        target_link_libraries(parameter_sweep mkl_rt)
        target_link_libraries(test_simulation mkl_rt)
        target_link_libraries(stroke_explore mkl_rt)
    endif()

endif()




find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(fic_explore Boost::program_options Boost::filesystem)
target_link_libraries(stroke_explore Boost::program_options Boost::filesystem)
target_link_libraries(parameter_sweep Boost::program_options)
target_link_libraries(test_simulation Boost::program_options)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tools)
include_directories(${CMAKE_SOURCE_DIR}/tools/optimlib)

find_package(Python3 3.8 COMPONENTS Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(fic_explore Python3::Python Python3::NumPy)
target_link_libraries(parameter_sweep Python3::Python Python3::NumPy)
target_link_libraries(test_simulation Python3::Python Python3::NumPy)
target_link_libraries(stroke_explore Python3::Python Python3::NumPy)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( fic_explore  ${ZLIB_LIBRARIES} )
target_link_libraries( parameter_sweep  ${ZLIB_LIBRARIES} )
target_link_libraries( test_simulation ${ZLIB_LIBRARIES} )
target_link_libraries( stroke_explore  ${ZLIB_LIBRARIES} )
