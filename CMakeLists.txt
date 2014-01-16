# CMakeLists.txt for the PBIL

cmake_minimum_required(VERSION 2.8)

# project settings
project(pbil)
set (pbil_VERSION_MAJOR 0)
set (pbil_VERSION_MINOR 1)


##### General Compilation Settings

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Add these standard paths to the search paths for FIND_LIBRARY
# to find libraries from these locations first
if(UNIX)
    set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH} /lib /usr/lib")
endif()

# --------------------------------------------------------------
# Indicate CMake 2.7 and above that we don't want to mix relative
#  and absolute paths in linker lib lists.
# Run "cmake --help-policy CMP0003" for more information.
# --------------------------------------------------------------
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif()

# ----------------------------------------------------------------------------
# Build static or dynamic libs?
# Default: dynamic libraries
# ----------------------------------------------------------------------------
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)")

# ----------------------------------------------------------------------------
# Include profiling?
# Default: NO
# ----------------------------------------------------------------------------
set(ENABLE_PROFILING OFF CACHE BOOL "Enable profiling in the GCC compiler (Add flags: -g -pg)")

# Profiling?
if(ENABLE_PROFILING)
    set(EXTRA_C_FLAGS_RELEASE "${EXTRA_C_FLAGS_RELEASE} -pg -g")
else()
    # Remove unreferenced functions: function level linking
    if(NOT APPLE)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -ffunction-sections")
    endif()
endif()

# LOGGING config
if(NOT LOG_ENGINE)
    set(LOG_ENGINE "BOOST_CUSTOM" CACHE STRING
        "choose a logging engine.  options are NO_LOG CERR BOOST_TRIVIAL BOOST_CUSTOM."
        FORCE)
endif(NOT LOG_ENGINE)
set(LOGGER_DEFINE "#define USE_LOGGER PBIL_LOGGING_${LOG_ENGINE}")


###### External Library Settings
#### Boost Library
#set(EXTRA_LIBS pthread)

set(Boost_USE_MULTITHREADED ON)
set(Boost_ADDITIONAL_VERSIONS "1.54" "1.54.0" "1.52" "1.52.0")
set(Boost_DETAILED_FAILURE_MSG ON)
find_path(BOOST_ROOT NAMES include/boost/circular_buffer.hpp
	HINTS /usr /usr/local /usr/local/share /opt/local ENV BOOST_HOME
	DOC "path to Boost root directory")
# must specify the components otherwise it's a header only usage
find_package( Boost REQUIRED COMPONENTS log thread system program_options)
if (Boost_FOUND)
    include_directories( ${Boost_INCLUDE_DIRS} )
    link_directories( ${Boost_LIBRARY_DIRS} )
    message(STATUS "Boost libs: " ${Boost_LIBRARIES})
endif(Boost_FOUND)
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})

# need this flag for boost log
add_definitions(-DBOOST_LOG_DYN_LINK)

#### MPI
set(USE_MPI ON CACHE BOOL "Build with MPI support")
if (USE_MPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
else(USE_MPI)
	#set(MPI_FOUND 0)
endif(USE_MPI)
if (MPI_FOUND)
    set(MPI_DEFINE "#define USE_MPI")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})
    
    message(STATUS "EXTRA libs: " ${EXTRA_LIBS})
else (MPI_FOUND)
    set(MPI_DEFINE "")    
endif (MPI_FOUND)

#### OpenMP
set(USE_OPENMP ON CACHE BOOL "Build with OpenMP support")
if (USE_OPENMP)
    include(FindOpenMP)

	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(USE_OPENMP)
    #set(OPENMP_FOUND 0)	
endif(USE_OPENMP)
if (OPENMP_FOUND)
    set(OPENMP_DEFINE "#define USE_OPENMP")
else (OPENMP_FOUND)
    set(OPENMP_DEFINE "")
endif (OPENMP_FOUND)


  

###### Executable and Libraries
# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Output directory for applications" )


#SET(LOGGER_DEFINE  CACHE STRING
#      "Choose logging engine, options are: N."
#      FORCE)
      
      
# Check if the user want to build test applications
option (BUILD_TEST_APPLICATIONS "Inform whether test applications should be built" OFF)
# Check if the user want to build sample applications
option (BUILD_EXAMPLE_APPLICATIONS "Inform whether sample applications should be built" OFF)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
    "${PROJECT_SOURCE_DIR}/config.hpp.in"
    "${PROJECT_BINARY_DIR}/config.hpp"
    @ONLY
    )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")


###### LIBRARY
#set(libsrcs CmdlineParser.cpp Communicator_I.cpp CommHandler_I.cpp Process.cpp PullCommHandler.cpp PushCommHandler.cpp RandomScheduler.cpp RoundRobinScheduler.cpp NullSinkAction.cpp DataBuffer.cpp MPIDataBuffer.cpp MPISendDataBuffer.cpp MPIRecvDataBuffer.cpp)
#set(libname runtime)
# Add rule to generate execution engine, as a shared library (.so)
#add_library("${libname}" SHARED ${libsrcs})
# Link opencv stuff to the execution engine

add_subdirectory(test)
#add_subdirectory(src/index)

# Install Settings


# Installer Settings


# Test Settings



