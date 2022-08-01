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
        include/simulator/monitors/raw.h
        include/simulator/process/bold_tvb.h
        src/simulator/process/bold_TVB.cpp
        include/simulator/process/bold_Stephan2007.h
        src/simulator/process/bold_Stephan2007.cpp
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
)

add_executable(parameter_sweep examples/parameter_sweep.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})
add_executable(test_simulation examples/test_simulation.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})
add_executable(plot_model examples/plot_model.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})
add_executable(plot_npz examples/plot_npz.cpp ${TVB_SOURCES})
add_executable(draft_cpp examples/Prepro_DecoEtAl2021_G_Optim.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES} include/tools/csv_tools.h src/tools/csv_tools.cpp examples/plot_model.cpp include/tools/threadpool.h)

find_package(MKL CONFIG)

if (MKL_FOUND)
    message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
    target_compile_options(draft_cpp PUBLIC $<TARGET_PROPERTY:MKL::mkl_core,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(draft_cpp INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    include_directories(${MKL_H})
    target_link_libraries(plot_model $<LINK_ONLY:MKL::mkl_core>)
    target_link_libraries(plot_model $<LINK_ONLY:MKL::mkl_intel_thread>)
endif()

if (MSVC)
    add_library(mkl_rt UNKNOWN IMPORTED)
    set_property(TARGET mkl_rt PROPERTY IMPORTED_LOCATION "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2022.0.3\\lib\\intel64\\mkl_rt.lib")
    target_link_libraries(parameter_sweep mkl_rt)
    target_link_libraries(test_simulation mkl_rt)
    target_link_libraries(plot_model mkl_rt)
endif()

if (UNIX)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g -ffast-math -DEIGEN_NO_DEBUG" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3  -DEIGEN_NO_DEBUG" CACHE STRING "" FORCE)
    # set(CMAKE_CXX_FLAGS_DEBUG  "-g -save-temps" CACHE STRING "" FORCE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(parameter_sweep Threads::Threads)

    if (MKL_FOUND)
        target_link_libraries(parameter_sweep mkl_rt)
        target_link_libraries(test_simulation mkl_rt)
        target_link_libraries(plot_model mkl_rt)
    endif()
endif()

find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

find_package(Boost 1.67.0 COMPONENTS program_options )
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(parameter_sweep Boost::program_options)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tools)

find_package(Python3 3.8 COMPONENTS Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(parameter_sweep Python3::Python Python3::NumPy)
target_link_libraries(test_simulation Python3::Python Python3::NumPy)
target_link_libraries(plot_model Python3::Python Python3::NumPy)
target_link_libraries(plot_npz Python3::Python Python3::NumPy)
target_link_libraries(draft_cpp Python3::Python Python3::NumPy)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( parameter_sweep  ${ZLIB_LIBRARIES} )
target_link_libraries( test_simulation ${ZLIB_LIBRARIES} )
target_link_libraries( plot_model ${ZLIB_LIBRARIES} )
target_link_libraries( plot_npz ${ZLIB_LIBRARIES} )
target_link_libraries( draft_cpp  ${ZLIB_LIBRARIES} )
