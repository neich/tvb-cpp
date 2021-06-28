cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELINTODEB  "-O2 -g" CACHE STRING "" FORCE)

project(tvb_cpp)
set(CMAKE_PROJECT_NAME "tvb-root-cpp")

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_VERBOSE_MAKEFILE ON)


set(TVB_SOURCES
        include/simulator/simulator.h
        src/simulator/simulator.cpp
        include/simulator/integrator.h
        include/simulator/model.h
        include/datatypes/connectivity.h
        include/definitions.h
        include/simulator/history.h
        src/simulator/history.cpp
        include/simulator/coupling.h
        src/simulator/coupling.cpp
        include/simulator/models/reduced_ww_ext.h
        src/simulator/models/reduced_ww_ext.cpp
        include/simulator/integrators/euler_deterministic.h
        src/simulator/integrators/euler_deterministic.cpp
        include/simulator/integrators/euler_stochastic.h
        src/simulator/integrators/euler_stochastic.cpp
        include/tools/cnpy.h
        src/tools/cnpy.cpp
        tools/fic/functions/balance_fic.h
        include/simulator/simulate.h
        src/simulator/simulate.cpp
        include/simulator/monitor.h
        include/simulator/monitors/raw.h
        include/simulator/monitors/bold_tvb.h
        src/simulator/monitors/bold_TVB.cpp
        include/simulator/monitors/bold_Stephan2007.h
        src/simulator/monitors/bold_Stephan2007.cpp
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

set_source_files_properties(src/simulator/models/reduced_ww_ext.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/integrators/euler_deterministic.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/integrators/euler_stochastic.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/coupling.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/history.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/simulate.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/simulator.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(src/simulator/noise.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/scipy/stats/stats.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/fic/functions/bold_filters.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/scipy/signal/signaltools.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/scipy/special/cephes/kolmogorov.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/numpy/numpy.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/scipy/linalg.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/numpy/core/numeric.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})
set_source_files_properties(tools/external/scipy/special/cephes/dd_real.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELINTODEB})

add_executable(test_simulation examples/test_simulation.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})
add_executable(plot_npz examples/plot_npz.cpp ${TVB_SOURCES})
add_executable(draft_cpp examples/Prepro_DecoEtAl2021_G_Optim.cpp ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})

find_package (Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tools)

find_package(Python3 COMPONENTS Development NumPy)
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(test_simulation Python3::Python Python3::NumPy)
target_link_libraries(plot_npz Python3::Python Python3::NumPy)
target_link_libraries(draft_cpp Python3::Python Python3::NumPy)

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( test_simulation ${ZLIB_LIBRARIES} )
target_link_libraries( plot_npz ${ZLIB_LIBRARIES} )
target_link_libraries( draft_cpp ${ZLIB_LIBRARIES} )

# find_package(MKL)
set(MKL_INCLUDE_DIR "/usr/include/mkl")
include_directories( ${MKL_INCLUDE_DIR} )

target_link_libraries( draft_cpp -lmkl_rt )
target_link_libraries( test_simulation -lmkl_rt )
