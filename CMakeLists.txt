cmake_minimum_required(VERSION 3.13)

project(tvb_cpp)
set(CMAKE_PROJECT_NAME "tvb-root-cpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (MSVC)
#    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "/Od" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELINTODEB  "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FS")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} /Ox")
endif()




set(CMAKE_CXX_FLAGS_RELINTODEB  "")

set(TVB_SOURCES
        include/simulator/simulator.h
        src/simulator/simulator.cpp
        include/simulator/integrator.h
        include/simulator/model.h
        include/datatypes/connectivity.h
        include/definitions.h
        include/simulator/coupling.h
        src/simulator/coupling.cpp
        include/simulator/models/reduced_ww_ext.h
        src/simulator/models/reduced_ww_ext.cpp
        include/simulator/models/zerlaut.h
        src/simulator/models/zerlaut.cpp
        include/simulator/models/montbrio.h
        src/simulator/models/montbrio.cpp
        include/simulator/integrators/euler_deterministic.h
        src/simulator/integrators/euler_deterministic.cpp
        include/simulator/integrators/euler_stochastic.h
        src/simulator/integrators/euler_stochastic.cpp
        include/tools/cnpy.h
        src/tools/cnpy.cpp
        include/tools/csv_tools.h
        src/tools/csv_tools.cpp
        tools/fic/functions/balance_fic.h
        include/simulator/simulate.h
        src/simulator/simulate.cpp
        include/simulator/monitor.h
        include/simulator/monitors/bold_tvb.h
        src/simulator/monitors/bold_tvb.cpp
        include/simulator/monitors/bold_Stephan2007.h
        src/simulator/monitors/bold_Stephan2007.cpp
        include/simulator/monitors/bold_BalloonWindkessel.h
        src/simulator/monitors/bold_BallonWindkessel.cpp
        include/datatypes/equation.h
        include/eigen/circshift.h
        include/tools/npz_tools.h
        src/tools/npz_tools.cpp
        )

set(TVB_TOOLS_SOURCES
        tools/fic/functions/balance_fic.cpp
        tools/fic/functions/observers/fc.h
        tools/fic/functions/observers/sw_fcd.h
        tools/fic/functions/observers/gbc.h
        tools/fic/functions/observers/sw_fcd.cpp
        tools/fic/functions/balance_fic.cpp
        tools/fic/functions/bold_filters.h
        tools/fic/functions/bold_filters.cpp
        tools/external/scipy/stats/stats.h
        tools/external/scipy/stats/stats.cpp
        tools/external/numpy/core/numeric.h
        tools/external/numpy/core/numeric.cpp
        tools/external/numpy/numpy.h
        tools/external/numpy/numpy.cpp
        tools/external/scipy/special/cephes/dd_real.cpp
        tools/external/scipy/special/cephes/kolmogorov.cpp
        tools/external/scipy/special/cephes/kolmogorov.h
        tools/external/scipy/stats/ksstats.h
        tools/external/scipy/stats/ksstats.cpp
        tools/external/scipy/signal/filter_design.h
        tools/external/scipy/signal/filter_design.cpp
        tools/external/scipy/signal/signaltools.h
        tools/external/scipy/signal/signaltools.cpp
        tools/external/scipy/linalg.h
        tools/external/scipy/linalg.cpp
        tools/external/eigen/eigen.h
        tools/fic/functions/observers/fc.h
        tools/fic/functions/g_optim.h
        tools/fic/functions/g_optim.cpp
        tools/fic/functions/simulate_fcd.h
        tools/fic/functions/simulate_fcd.cpp
        tools/load_or_compute.h
        include/tools/npy.h)

add_library(tvb_cpp SHARED ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})

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
    target_link_libraries(tvb_cpp Threads::Threads)

    if (MKL_FOUND)
        target_link_libraries(tvb_cpp  mkl_rt)
    endif()
    if (BLAS_FOUND)
        target_link_libraries(tvb_cpp mkl_rt)
    endif()

endif()


find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(tvb_cpp Boost::program_options Boost::filesystem)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tools)
include_directories(${CMAKE_SOURCE_DIR}/tools/optimlib)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( tvb_cpp  ${ZLIB_LIBRARIES} )
