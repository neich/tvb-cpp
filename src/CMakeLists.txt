cmake_minimum_required(VERSION 3.13)

project(tvb-cpp)
set(CMAKE_PROJECT_NAME "tvb-cpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (MSVC)
#    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "/Od" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELINTODEB  "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FS")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} /Ox")
endif()


set(TVB_SOURCES
        tvb-cpp/simulator/simulator.h
        tvb-cpp/simulator/simulator.cpp
        tvb-cpp/simulator/integrator.h
        tvb-cpp/simulator/model.h
        tvb-cpp/datatypes/connectivity.h
        tvb-cpp/definitions.h
        tvb-cpp/simulator/coupling.h
        tvb-cpp/simulator/coupling.cpp
        tvb-cpp/simulator/models/reduced_ww_ext.h
        tvb-cpp/simulator/models/reduced_ww_ext.cpp
        tvb-cpp/simulator/models/zerlaut.h
        tvb-cpp/simulator/models/zerlaut.cpp
        tvb-cpp/simulator/models/montbrio.h
        tvb-cpp/simulator/models/montbrio.cpp
        tvb-cpp/simulator/integrators/euler_deterministic.h
        tvb-cpp/simulator/integrators/euler_deterministic.cpp
        tvb-cpp/simulator/integrators/euler_stochastic.h
        tvb-cpp/simulator/integrators/euler_stochastic.cpp
        tvb-cpp/tools/cnpy.h
        tvb-cpp/tools/cnpy.cpp
        tvb-cpp/tools/csv_tools.h
        tvb-cpp/tools/csv_tools.cpp
        tvb-cpp/tools/algo/fic/functions/balance_fic.h
        tvb-cpp/simulator/simulate.h
        tvb-cpp/simulator/simulate.cpp
        tvb-cpp/simulator/monitor.h
        tvb-cpp/simulator/monitors/bold_tvb.h
        tvb-cpp/simulator/monitors/bold_tvb.cpp
        tvb-cpp/simulator/monitors/bold_Stephan2007.h
        tvb-cpp/simulator/monitors/bold_Stephan2007.cpp
        tvb-cpp/simulator/monitors/bold_BalloonWindkessel.h
        tvb-cpp/simulator/monitors/bold_BallonWindkessel.cpp
        tvb-cpp/datatypes/equation.h
        tvb-cpp/eigen/circshift.h
        tvb-cpp/tools/npz_tools.h
        tvb-cpp/tools/npz_tools.cpp
        )

set(TVB_TOOLS_SOURCES
        tvb-cpp/tools/algo/fic/functions/balance_fic.cpp
        tvb-cpp/tools/algo/fic/functions/observers/fc.h
        tvb-cpp/tools/algo/fic/functions/observers/sw_fcd.h
        tvb-cpp/tools/algo/fic/functions/observers/gbc.h
        tvb-cpp/tools/algo/fic/functions/observers/sw_fcd.cpp
        tvb-cpp/tools/algo/fic/functions/balance_fic.cpp
        tvb-cpp/tools/algo/fic/functions/bold_filters.h
        tvb-cpp/tools/algo/fic/functions/bold_filters.cpp
        tvb-cpp/tools/algo/external/scipy/stats/stats.h
        tvb-cpp/tools/algo/external/scipy/stats/stats.cpp
        tvb-cpp/tools/algo/external/numpy/core/numeric.h
        tvb-cpp/tools/algo/external/numpy/core/numeric.cpp
        tvb-cpp/tools/algo/external/numpy/numpy.h
        tvb-cpp/tools/algo/external/numpy/numpy.cpp
        tvb-cpp/tools/algo/external/scipy/special/cephes/dd_real.cpp
        tvb-cpp/tools/algo/external/scipy/special/cephes/kolmogorov.cpp
        tvb-cpp/tools/algo/external/scipy/special/cephes/kolmogorov.h
        tvb-cpp/tools/algo/external/scipy/stats/ksstats.h
        tvb-cpp/tools/algo/external/scipy/stats/ksstats.cpp
        tvb-cpp/tools/algo/external/scipy/signal/filter_design.h
        tvb-cpp/tools/algo/external/scipy/signal/filter_design.cpp
        tvb-cpp/tools/algo/external/scipy/signal/signaltools.h
        tvb-cpp/tools/algo/external/scipy/signal/signaltools.cpp
        tvb-cpp/tools/algo/external/scipy/linalg.h
        tvb-cpp/tools/algo/external/scipy/linalg.cpp
        tvb-cpp/tools/algo/external/eigen/eigen.h
        tvb-cpp/tools/algo/fic/functions/observers/fc.h
        tvb-cpp/tools/algo/fic/functions/g_optim.h
        tvb-cpp/tools/algo/fic/functions/g_optim.cpp
        tvb-cpp/tools/algo/fic/functions/simulate_fcd.h
        tvb-cpp/tools/algo/fic/functions/simulate_fcd.cpp
        tvb-cpp/tools/algo/load_or_compute.h
        tvb-cpp/tools/npy.h
        tvb-cpp/tools/json.h)

add_library(${CMAKE_PROJECT_NAME} SHARED ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})

find_package(MKL CONFIG)

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELDEBPROFILE  "-O2 -march=native -pg -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -march=native -g -DEIGEN_NO_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DEIGEN_NO_DEBUG -DNDEBUG")
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

