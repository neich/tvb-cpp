cmake_minimum_required(VERSION 3.13)

project(tvb-root-cpp)
set(CMAKE_PROJECT_NAME "tvb-root-cpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (MSVC)
#    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "/Od" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELINTODEB  "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FS")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} /Ox")
endif()


set(TVB_SOURCES
        tvb-root-cpp/simulator/simulator.h
        tvb-root-cpp/simulator/simulator.cpp
        tvb-root-cpp/simulator/integrator.h
        tvb-root-cpp/simulator/model.h
        tvb-root-cpp/datatypes/connectivity.h
        tvb-root-cpp/definitions.h
        tvb-root-cpp/simulator/coupling.h
        tvb-root-cpp/simulator/coupling.cpp
        tvb-root-cpp/simulator/models/reduced_ww_ext.h
        tvb-root-cpp/simulator/models/reduced_ww_ext.cpp
        tvb-root-cpp/simulator/models/zerlaut.h
        tvb-root-cpp/simulator/models/zerlaut.cpp
        tvb-root-cpp/simulator/models/montbrio.h
        tvb-root-cpp/simulator/models/montbrio.cpp
        tvb-root-cpp/simulator/integrators/euler_deterministic.h
        tvb-root-cpp/simulator/integrators/euler_deterministic.cpp
        tvb-root-cpp/simulator/integrators/euler_stochastic.h
        tvb-root-cpp/simulator/integrators/euler_stochastic.cpp
        tvb-root-cpp/tools/cnpy.h
        tvb-root-cpp/tools/cnpy.cpp
        tvb-root-cpp/tools/csv_tools.h
        tvb-root-cpp/tools/csv_tools.cpp
        tvb-root-cpp/tools/algo/fic/functions/balance_fic.h
        tvb-root-cpp/simulator/simulate.h
        tvb-root-cpp/simulator/simulate.cpp
        tvb-root-cpp/simulator/monitor.h
        tvb-root-cpp/simulator/monitors/bold_tvb.h
        tvb-root-cpp/simulator/monitors/bold_tvb.cpp
        tvb-root-cpp/simulator/monitors/bold_Stephan2007.h
        tvb-root-cpp/simulator/monitors/bold_Stephan2007.cpp
        tvb-root-cpp/simulator/monitors/bold_BalloonWindkessel.h
        tvb-root-cpp/simulator/monitors/bold_BallonWindkessel.cpp
        tvb-root-cpp/datatypes/equation.h
        tvb-root-cpp/eigen/circshift.h
        tvb-root-cpp/tools/npz_tools.h
        tvb-root-cpp/tools/npz_tools.cpp
        )

set(TVB_TOOLS_SOURCES
        tvb-root-cpp/tools/algo/fic/functions/balance_fic.cpp
        tvb-root-cpp/tools/algo/fic/functions/observers/fc.h
        tvb-root-cpp/tools/algo/fic/functions/observers/sw_fcd.h
        tvb-root-cpp/tools/algo/fic/functions/observers/gbc.h
        tvb-root-cpp/tools/algo/fic/functions/observers/sw_fcd.cpp
        tvb-root-cpp/tools/algo/fic/functions/balance_fic.cpp
        tvb-root-cpp/tools/algo/fic/functions/bold_filters.h
        tvb-root-cpp/tools/algo/fic/functions/bold_filters.cpp
        tvb-root-cpp/tools/algo/external/scipy/stats/stats.h
        tvb-root-cpp/tools/algo/external/scipy/stats/stats.cpp
        tvb-root-cpp/tools/algo/external/numpy/core/numeric.h
        tvb-root-cpp/tools/algo/external/numpy/core/numeric.cpp
        tvb-root-cpp/tools/algo/external/numpy/numpy.h
        tvb-root-cpp/tools/algo/external/numpy/numpy.cpp
        tvb-root-cpp/tools/algo/external/scipy/special/cephes/dd_real.cpp
        tvb-root-cpp/tools/algo/external/scipy/special/cephes/kolmogorov.cpp
        tvb-root-cpp/tools/algo/external/scipy/special/cephes/kolmogorov.h
        tvb-root-cpp/tools/algo/external/scipy/stats/ksstats.h
        tvb-root-cpp/tools/algo/external/scipy/stats/ksstats.cpp
        tvb-root-cpp/tools/algo/external/scipy/signal/filter_design.h
        tvb-root-cpp/tools/algo/external/scipy/signal/filter_design.cpp
        tvb-root-cpp/tools/algo/external/scipy/signal/signaltools.h
        tvb-root-cpp/tools/algo/external/scipy/signal/signaltools.cpp
        tvb-root-cpp/tools/algo/external/scipy/linalg.h
        tvb-root-cpp/tools/algo/external/scipy/linalg.cpp
        tvb-root-cpp/tools/algo/external/eigen/eigen.h
        tvb-root-cpp/tools/algo/fic/functions/observers/fc.h
        tvb-root-cpp/tools/algo/fic/functions/g_optim.h
        tvb-root-cpp/tools/algo/fic/functions/g_optim.cpp
        tvb-root-cpp/tools/algo/fic/functions/simulate_fcd.h
        tvb-root-cpp/tools/algo/fic/functions/simulate_fcd.cpp
        tvb-root-cpp/tools/algo/load_or_compute.h
        tvb-root-cpp/tools/npy.h)

add_library(${CMAKE_PROJECT_NAME} SHARED ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})

find_package(MKL CONFIG)

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -march=native -g -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O2 -march=native -DEIGEN_NO_DEBUG -DNDEBUG")
    # set(CMAKE_CXX_FLAGS_DEBUG  "-g -save-temps" CACHE STRING "" FORCE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(${CMAKE_PROJECT_NAME} Threads::Threads)

    if (MKL_FOUND)
        add_definitions( -DEIGEN_USE_MKL_ALL )
        include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(${CMAKE_PROJECT_NAME}  $<LINK_ONLY:MKL::MKL>)
    endif()
endif()


find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} Boost::program_options Boost::filesystem)

include_directories(${PROJECT_SOURCE_DIR})

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( ${CMAKE_PROJECT_NAME}  ${ZLIB_LIBRARIES} )

