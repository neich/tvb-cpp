cmake_minimum_required(VERSION 3.15...3.22)

project(tvb_root_cpp_python VERSION "0.0.6")
set(CMAKE_CXX_STANDARD 20)

#if(SKBUILD)
#  # Scikit-Build does not add your site-packages to the search path
#  # automatically, so we need to add it _or_ the pybind11 specific directory
#  # here.
#  execute_process(
#    COMMAND "${PYTHON_EXECUTABLE}" -c
#            "import pybind11; print(pybind11.get_cmake_dir())"
#    OUTPUT_VARIABLE _tmp_dir
#    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
#  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
#endif()

# Now we can find pybind11
add_subdirectory(pybind11)

if (MSVC)
  #    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "/Od" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:preprocessor /std:c++20")
endif()

set(TVB_SOURCES
        ../src/tvb-cpp/simulator/simulator.h
        ../src/tvb-cpp/simulator/simulator.cpp
        ../src/tvb-cpp/simulator/integrator.h
        ../src/tvb-cpp/simulator/model.h
        ../src/tvb-cpp/datatypes/connectivity.h
        ../src/tvb-cpp/definitions.h
        ../src/tvb-cpp/simulator/coupling.h
        ../src/tvb-cpp/simulator/coupling.cpp
        ../src/tvb-cpp/simulator/models/reduced_ww_ext.h
        ../src/tvb-cpp/simulator/models/reduced_ww_ext.cpp
        ../src/tvb-cpp/simulator/models/zerlaut.h
        ../src/tvb-cpp/simulator/models/zerlaut.cpp
        ../src/tvb-cpp/simulator/models/montbrio.h
        ../src/tvb-cpp/simulator/models/montbrio.cpp
        ../src/tvb-cpp/simulator/integrators/euler_deterministic.h
        ../src/tvb-cpp/simulator/integrators/euler_deterministic.cpp
        ../src/tvb-cpp/simulator/integrators/euler_stochastic.h
        ../src/tvb-cpp/simulator/integrators/euler_stochastic.cpp
        ../src/tvb-cpp/tools/algo/fic/functions/balance_fic.h
        ../src/tvb-cpp/simulator/simulate.h
        ../src/tvb-cpp/simulator/simulate.cpp
        ../src/tvb-cpp/simulator/monitor.h
        ../src/tvb-cpp/simulator/bold/bold_tvb.h
        ../src/tvb-cpp/simulator/bold/bold_tvb.cpp
        ../src/tvb-cpp/simulator/bold/bold_Stephan2007.h
        ../src/tvb-cpp/simulator/bold/bold_Stephan2007.cpp
        ../src/tvb-cpp/simulator/bold/bold_BalloonWindkessel.h
        ../src/tvb-cpp/simulator/bold/bold_BallonWindkessel.cpp
        ../src/tvb-cpp/datatypes/equation.h
        ../src/tvb-cpp/eigen/circshift.h
        )

set(TVB_TOOLS_SOURCES
        ../src/tvb-cpp/tools/algo/fic/functions/balance_fic.cpp
        ../src/tvb-cpp/tools/observers/fc.h
        ../src/tvb-cpp/tools/observers/sw_fcd.h
        ../src/tvb-cpp/tools/observers/gbc.h
        ../src/tvb-cpp/tools/observers/sw_fcd.cpp
        ../src/tvb-cpp/tools/algo/fic/functions/balance_fic.cpp
        ../src/tvb-cpp/tools/bold_filters.h
        ../src/tvb-cpp/tools/bold_filters.cpp
        ../src/tvb-cpp/tools/scipy/stats/stats.h
        ../src/tvb-cpp/tools/scipy/stats/stats.cpp
        ../src/tvb-cpp/tools/numpy/core/numeric.h
        ../src/tvb-cpp/tools/numpy/core/numeric.cpp
        ../src/tvb-cpp/tools/numpy/numpy.h
        ../src/tvb-cpp/tools/numpy/numpy.cpp
        ../src/tvb-cpp/tools/scipy/special/cephes/dd_real.cpp
        ../src/tvb-cpp/tools/scipy/special/cephes/kolmogorov.cpp
        ../src/tvb-cpp/tools/scipy/special/cephes/kolmogorov.h
        ../src/tvb-cpp/tools/scipy/stats/ksstats.h
        ../src/tvb-cpp/tools/scipy/stats/ksstats.cpp
        ../src/tvb-cpp/tools/scipy/signal/filter_design.h
        ../src/tvb-cpp/tools/scipy/signal/filter_design.cpp
        ../src/tvb-cpp/tools/scipy/signal/signaltools.h
        ../src/tvb-cpp/tools/scipy/signal/signaltools.cpp
        ../src/tvb-cpp/tools/scipy/linalg.h
        ../src/tvb-cpp/tools/scipy/linalg.cpp
        ../src/tvb-cpp/tools/scipy/fft.h
        ../src/tvb-cpp/tools/scipy/fft.cpp
        ../src/tvb-cpp/tools/eigen/eigen.h
        ../src/tvb-cpp/tools/observers/fc.h
        ../src/tvb-cpp/tools/algo/load_or_compute.h)



find_package(Eigen3 3.4 REQUIRED)
include_directories($<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>)


find_package(Boost 1.67.0 COMPONENTS program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})

find_package( ZLIB )

include_directories(${ZLIB_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/../src)

pybind11_add_module(_core MODULE src/main.cpp src/api.cpp src/api.h ${TVB_SOURCES} ${TVB_TOOLS_SOURCES})

if (UNIX)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g -DEIGEN_NO_DEBUG" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -ffast-math -DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
  # set(CMAKE_CXX_FLAGS_DEBUG  "-g -save-temps" CACHE STRING "" FORCE)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(_core PRIVATE Threads::Threads)
endif()

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_core PRIVATE Boost::program_options Boost::filesystem)
target_link_libraries(_core PRIVATE ${ZLIB_LIBRARIES})

install(TARGETS _core DESTINATION .)
