# CMakeLists.txt for the BLISS

cmake_minimum_required(VERSION 2.6)

# project settings
project(bliss)
set (bliss_VERSION_MAJOR 0)
set (bliss_VERSION_MINOR 1)

INCLUDE(CMakeDependentOption)

##### General Compilation Settings

### REQUIRE C++11
# Initialize CXXFLAGS.
add_definitions(-Wall -std=c++11 -fPIC)


### tool to speed up compilation
if(CMAKE_VERSION VERSION_GREATER 2.8.3)
    set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/cotire/CMake")
    include(cotire)
else()
    function(cotire module_names)
    endfunction(cotire)  
endif()



# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
#    message(FATAL_ERROR "${PROJECT_NAME} currently does not support Intel compiler due to constexpr not allowing math function calls")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE ICPC_VERSION)
  if (NOT (ICPC_VERSION VERSION_GREATER 13.0 OR ICPC_VERSION VERSION_EQUAL 13.0))
    message(FATAL_ERROR "${PROJECT_NAME} requires icpc 13.0 or greater.")
  endif ()
  SET(SUPPORTS_SANITIZER OFF CACHE INTERNAL "compiler supports sanitizer?")
  SET(SUPPORTS_COVERAGE OFF CACHE INTERNAL "compiler supports coverage?")

  SET(EXTRA_WARNING_FLAGS "" CACHE INTERNAL "extra compiler warning flags")
  SET(TYPE_CONVERSION_WARNING_FLAGS "" CACHE INTERNAL "compiler flags to check type conversions")
  SET(SUGGESTION_WARNING_FLAGS CACHE INTERNAL "compiler flags that suggest keywords for better type resolutions")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
  if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
    message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
  endif ()
  SET(SUPPORTS_SANITIZER ON CACHE INTERNAL "compiler supports sanitizer?")
  SET(SUPPORTS_COVERAGE ON CACHE INTERNAL "compiler supports coverage?")
  
# -Wzero-as-null-pointer-constant : causes a lot of errors in system headers.
  SET(EXTRA_WARNING_FLAGS "-Wextra -Wuninitialized -Wswitch-default -Wundef")
  SET(TYPE_CONVERSION_WARNING_FLAGS "-Wdouble-promotion -Wconversion -Wsign-conversion -Wcast-qual -Wuseless-cast -Wfloat-equal")
  SET(SUGGESTION_WARNING_FLAGS "-Wsuggest-override -Wsuggest-final-types -Wsuggest-final-methods")
  
  
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "${PROJECT_NAME} currently does not support clang due to include path issues")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
#  set(_CMAKE_TOOLCHAIN_PREFIX "llvm-")
#  include_directories("/usr/local/lib/clang/3.4/include")
  SET(SUPPORTS_SANITIZER ON CACHE INTERNAL "compiler supports sanitizer?")
  SET(SUPPORTS_COVERAGE ON CACHE INTERNAL "compiler supports coverage?")

  SET(EXTRA_WARNING_FLAGS "" CACHE INTERNAL "extra compiler warning flags")
  SET(TYPE_CONVERSION_WARNING_FLAGS "" CACHE INTERNAL "compiler flags to check type conversions")
  SET(SUGGESTION_WARNING_FLAGS CACHE INTERNAL "compiler flags that suggest keywords for better type resolutions")
  
else ()
  message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()
### from http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake


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

  if(CMAKE_VERSION VERSION_GREATER 2.6.4)
    cmake_policy(SET CMP0012 NEW)
  endif()
endif()

# ----------------------------------------------------------------------------
# Build static or dynamic libs?
# Default: dynamic libraries
# ----------------------------------------------------------------------------
OPTION(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" ON)

# ----------------------------------------------------------------------------
# Include profiling?
# Default: NO
# ----------------------------------------------------------------------------
OPTION(ENABLE_PROFILING "Enable profiling in the GCC compiler (Add flags: -g -pg)" OFF)

# Profiling?
if(ENABLE_PROFILING)
  add_definitions(-pg -g)
else()
  # Remove unreferenced functions: function level linking
  if(NOT APPLE)
    add_definitions(-ffunction-sections)
  endif()
endif()

OPTION(ENABLE_EXTRA_COMPILER_WARNINGS "Enable compiler to generate extra warnings?" OFF)
OPTION(ENABLE_TYPE_CONVERSION_WARNINGS "Enable compiler to generate extra warnings for type conversions?" OFF)
OPTION(ENABLE_SUGGESTION_WARNINGS "Enable compiler to generate suggestions for polymorphic type analysis?" OFF)


###### LOGGING ENGINE configuration
SET(LOGGER_VERBOSITY "WARNING" CACHE STRING "Logger verbosity.  Valid entries are ERROR, WARNING, INFO, DEBUG, TRACE.  Each level includes all levels with smaller value")  
#function(update_logger_verbosity varname cacheaccess varval) 
#    SET(LOGGER_VERBOSITY_DEFINE "#define LOGGER_VERBOSITY BLISS_LOGGER_VERBOSITY_${LOGGER_VERBOSITY}")
#endfunction(update_logger_verbosity varname cacheaccess varval)
#variable_watch(LOGGER_VERBOSITY update_logger_verbosity)

function(update_logging_engine varname cacheaccess varval) 
#       message(STATUS "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")


    if (${cacheaccess} STREQUAL READ_ACCESS AND
        ${varname} STREQUAL USE_OPENMP)

            
        if (USE_OPENMP)

            # restrict log engine to no_log or printf
            if (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "PRINTF")
                
                set(LOG_ENGINE ${LOG_ENGINE} CACHE STRING
                "choose a logging engine.  options are NO_LOG PRINTF." FORCE)
                
            else (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "PRINTF")

#                message(STATUS "OMP ENABLED.  Default Log Engine set to NO_LOG")
                
                set(LOG_ENGINE "NO_LOG" CACHE STRING
                "choose a logging engine.  options are NO_LOG PRINTF." FORCE)
                
            endif (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "PRINTF")
            
        else (USE_OPENMP)
            # OMP debugging is not on.  so log engine choice depends on whether boost logging is enabled.
        
          if (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "CERR" OR
                LOG_ENGINE STREQUAL "PRINTF" OR
                LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")
                
                
            set(LOG_ENGINE ${LOG_ENGINE} CACHE STRING
            "choose a logging engine.  options are NO_LOG PRINTF CERR BOOST_TRIVIAL BOOST_CUSTOM." FORCE)
          else (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "CERR" OR
                LOG_ENGINE STREQUAL "PRINTF" OR
                LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")
                
#                message(STATUS "OMP DISABLED.  Default Log Engine set to NO_LOG")
                
            set(LOG_ENGINE "NO_LOG" CACHE STRING
            "choose a logging engine.  options are NO_LOG PRINTF CERR BOOST_TRIVIAL BOOST_CUSTOM." FORCE)
          endif (LOG_ENGINE STREQUAL "NO_LOG" OR
                LOG_ENGINE STREQUAL "CERR" OR
                LOG_ENGINE STREQUAL "PRINTF" OR
                LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")

        endif(USE_OPENMP)

    else(${cacheaccess} STREQUAL READ_ACCESS AND




    
        ${varname} STREQUAL USE_OPENMP)
        
       message(WARNING "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")
        
    endif(${cacheaccess} STREQUAL READ_ACCESS AND
        ${varname} STREQUAL USE_OPENMP)
        
        
    # set the log engine in the logger_config.hpp file.
    set(LOGGER_DEFINE "#define USE_LOGGER BLISS_LOGGING_${LOG_ENGINE}")
#    configure_file(
#        "${PROJECT_SOURCE_DIR}/src/config/logger_config.hpp.in"
#        "${PROJECT_BINARY_DIR}/logger_config.hpp"
#        @ONLY
#        )
#        
endfunction(update_logging_engine varname cacheaccess varval)



# if OMP_DEBUGGING is enabled, default(none) is set in the omp pragmas.   this means that only NO_LOG and PRINTF logging are compatible.


CMAKE_DEPENDENT_OPTION(OPENMP_STRICT_SCOPING "Enable OpenMP debugging (turns on default(none) in omp pragma)" ON
                        "OMP_DEBUGGING" OFF)
mark_as_advanced(OPENMP_STRICT_SCOPING)
if (OPENMP_STRICT_SCOPING)
    add_definitions(-DOMP_DEBUG)
endif(OPENMP_STRICT_SCOPING)
    




###### External Library Settings
#### Boost Library
#set(EXTRA_LIBS pthread)
if (LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_ADDITIONAL_VERSIONS "1.54" "1.54.0")
    set(Boost_DETAILED_FAILURE_MSG ON)
    find_path(BOOST_ROOT NAMES include/boost/circular_buffer.hpp
      HINTS /usr /usr/local /usr/local/share /opt/local ENV BOOST_HOME
      DOC "path to Boost root directory")
    # must specify the components otherwise it's a header only usage
    find_package(Boost REQUIRED COMPONENTS system program_options log)
    if (Boost_FOUND)
      include_directories( ${Boost_INCLUDE_DIRS} )
      link_directories( ${Boost_LIBRARY_DIRS} )
      set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})
      add_definitions(-DBOOST_LOG_DYN_LINK)
    else(Boost_FOUND)
        unset(Boost_USE_MULTITHREADED)
        unset(Boost_ADDITIONAL_VERSIONS)
        unset(Boost_DETAILED_FAILURE_MSG)
        unset(BOOST_ROOT CACHE)
        
        set(LOG_ENGINE "NO_LOG" CACHE STRING
            "choose a logging engine.  options are NO_LOG PRINTF CERR BOOST_TRIVIAL BOOST_CUSTOM." FORCE)
        message(WARNING "Did not find boost.  Default Log Engine set to NO_LOG")
        
        set(LOGGER_DEFINE "#define USE_LOGGER BLISS_LOGGING_${LOG_ENGINE}")
#        configure_file(
#            "${PROJECT_SOURCE_DIR}/src/config/logger_config.hpp.in"
#            "${PROJECT_BINARY_DIR}/logger_config.hpp"
#            @ONLY
#            )
    endif(Boost_FOUND)
    
    # need this flag for boost log
else(LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")
    unset(Boost_USE_MULTITHREADED)
    unset(Boost_ADDITIONAL_VERSIONS)
    unset(Boost_DETAILED_FAILURE_MSG)
    unset(Boost_FOUND)
    unset(BOOST_ROOT CACHE)
    unset(Boost_DIR CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_LIBRARY_DIR CACHE)
    unset(Boost_LOG_LIBRARY_DEBUG CACHE)
    unset(Boost_LOG_LIBRARY_RELEASE CACHE)
    unset(Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG CACHE)
    unset(Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE CACHE)
    unset(Boost_SYSTEM_LIBRARY_DEBUG CACHE)
    unset(Boost_SYSTEM_LIBRARY_RELEASE CACHE)

endif(LOG_ENGINE STREQUAL "BOOST_TRIVIAL" OR
                LOG_ENGINE STREQUAL "BOOST_CUSTOM")


#### MPI
OPTION(USE_MPI "Build with MPI support" ON)
if (USE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  message(STATUS "Found MPI:")
  message(STATUS "    headers: ${MPI_INCLUDE_PATH}")
  message(STATUS "    libs:    ${MPI_LIBRARIES}")
  add_definitions(${MPI_COMPILE_FLAGS})
  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_EXE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
  set(CMAKE_MODULE_LINK_FLAGS "${CMAKE_MODULE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
  set(CMAKE_SHARED_LINK_FLAGS "${CMAKE_SHARED_LINK_FLAGS} ${MPI_LINK_FLAGS}")
  set(CMAKE_STATIC_LINK_FLAGS "${CMAKE_STATIC_LINK_FLAGS} ${MPI_LINK_FLAGS}")
else(USE_MPI)
  set(MPI_FOUND 0)
endif(USE_MPI)

if (MPI_FOUND)
  set(MPI_DEFINE "#define USE_MPI")
  set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})
  set(ENV{OMPI_CC} ${CMAKE_C_COMPILER})
  set(ENV{OMPI_CXX} ${CMAKE_CXX_COMPILER})
  mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY)
  
else (MPI_FOUND)
  set(MPI_DEFINE "")
  message(WARNING "Not using MPI")
endif (MPI_FOUND)

#### OpenMP
include(FindOpenMP)

if (OPENMP_FOUND)
    message(STATUS "Found OpenMP")
endif(OPENMP_FOUND)

CMAKE_DEPENDENT_OPTION(USE_OPENMP "Build with OpenMP support" ON
                        "OPENMP_FOUND" OFF)
variable_watch(USE_OPENMP update_logging_engine)

if (USE_OPENMP)
    # add OpenMP flags to compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(OPENMP_DEFINE "#define USE_OPENMP")
else (USE_OPENMP)
    set(OPENMP_DEFINE "")
endif (USE_OPENMP)

    


###### Doxygen documentation

find_package(Doxygen)

if(DOXYGEN_FOUND)
  message(STATUS "Found Doxygen")
  message(STATUS "    use `make doc` to build documentation")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)



###### TESTING

OPTION(ENABLE_TESTING "Enable different testing mechanisms including coverage, thread and address sanitizers, and google test" OFF)
CMAKE_DEPENDENT_OPTION(ENABLE_COVERAGE "Enable Code Coverage Reporting" OFF
                        "ENABLE_TESTING;SUPPORTS_COVERAGE" OFF)
CMAKE_DEPENDENT_OPTION(ENABLE_SANITIZER "Enable Sanitizer.  Choose via SANITIZER_STYLE" OFF
                        "ENABLE_TESTING;SUPPORTS_SANITIZER" OFF)
                        
include("${CMAKE_SOURCE_DIR}/cmake/GCCFilter.cmake")

## LOGGING ENGINE configuration
#function(update_test_options varname cacheaccess varval) 
##       message(STATUS "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")
#
#
#    if (${cacheaccess} STREQUAL READ_ACCESS AND
#        ${varname} STREQUAL ENABLE_TESTING)
#
#            
#        if (ENABLE_TESTING)
#            if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#              set(ENABLE_COVERAGE ${ENABLE_COVERAGE} CACHE BOOL "Enable Code Coverage reporting" FORCE)
#              set(ENABLE_SANITIZER ${ENABLE_SANITIZER} CACHE BOOL "Enable Sanitizer.  Choose via SANITIZER_STYLE" FORCE)              
#            elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#              set(ENABLE_COVERAGE ${ENABLE_COVERAGE} CACHE BOOL "Enable Code Coverage reporting" FORCE)
#              set(ENABLE_SANITIZER ${ENABLE_SANITIZER} CACHE BOOL "Enable Sanitizer.  Choose via SANITIZER_STYLE" FORCE)
#            else()
#              set(ENABLE_COVERAGE OFF CACHE INTERNAL "DISABLE Code Coverage Reporting because we are using Intel Compiler - don't know profiling options." FORCE) 
#              set(ENABLE_SANITIZER OFF CACHE INTERNAL "Disabled Sanitizer as we are using Intel compiler." FORCE)
#            endif()
#            
#        else()
#              set(ENABLE_COVERAGE OFF CACHE INTERNAL "DISABLE Code Coverage Reporting because we are using Intel Compiler - don't know profiling options." FORCE) 
#              set(ENABLE_SANITIZER OFF CACHE INTERNAL "Disabled Sanitizer as testing is turned off." FORCE)
#        endif()
#
#        if (ENABLE_SANITIZER)
#            set(SANITIZER_STYLE ${SANITIZER_STYLE} CACHE STRING "Any Compiler supported Sanitizer style.  e.g. thread, address" FORCE)
#        else()
#            set(SANITIZER_STYLE "" CACHE INTERNAL "Disabled Sanitizer as we are NOT enabling sanitizer." FORCE)
#        endif()
#            
#                
#
#
#    else(${cacheaccess} STREQUAL READ_ACCESS AND
#        ${varname} STREQUAL ENABLE_TESTING)
#        
#       message(WARNING "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")
#        
#    endif(${cacheaccess} STREQUAL READ_ACCESS AND
#        ${varname} STREQUAL ENABLE_TESTING)
#                
#endfunction(update_test_options varname cacheaccess varval)
## Testing can be turned off: (for now the default is: ON)
#variable_watch(ENABLE_TESTING update_test_options)

# LOGGING ENGINE configuration
function(update_sanitizer_options varname cacheaccess varval) 
#       message(STATUS "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")


    if (${cacheaccess} STREQUAL READ_ACCESS AND
        ${varname} STREQUAL ENABLE_SANITIZER)

            
        if (${varval})
            set(SANITIZER_STYLE ${SANITIZER_STYLE} CACHE STRING "Any Compiler supported Sanitizer style.  e.g. thread, address" FORCE)
        else(${varval})
            set(SANITIZER_STYLE "" CACHE INTERNAL "Disabled Sanitizer as we are NOT enabling sanitizer." FORCE)
        endif(${varval})
            
    else(${cacheaccess} STREQUAL READ_ACCESS AND
        ${varname} STREQUAL ENABLE_SANITIZER)
        
       message(WARNING "${varname} changed to ${varval} with unhandled access: ${cacheaccess}.")
        
    endif(${cacheaccess} STREQUAL READ_ACCESS AND
        ${varname} STREQUAL ENABLE_SANITIZER)
                
endfunction(update_sanitizer_options varname cacheaccess varval)


variable_watch(ENABLE_SANITIZER update_sanitizer_options)



###### Google test library
# require "git submodule init" and "git submodule update"
if (ENABLE_TESTING)

    if (ENABLE_COVERAGE)
      # set flags for coverage test
      message(STATUS "Code coverage reporting enabled")
      add_definitions(-fprofile-arcs -ftest-coverage -g)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
      set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -g")
    endif()



  message(STATUS "Testing enabled")
  # include/load the Google Test framework
  set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
#  set_property(DIRECTORY ${EXT_PROJECTS_DIR}/gtest APPEND PROPERTY CMAKE_CXX_FLAGS -fPIC)
#  set_property(DIRECTORY ${EXT_PROJECTS_DIR}/gtest APPEND PROPERTY CMAKE_C_FLAGS -fPIC)
  add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
  include_directories(${EXT_PROJECTS_DIR}/gtest/include)

  # enable the CMAKE testing framework
  enable_testing()

  # Set paths outside of the function so the function can
  # be called from anywhere in the source tree
  #   set path for output of testing binaries
  set(TEST_BINARY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/test)
  #   set path for output of test reports (XML format)
  set(TEST_XML_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Testing)


  function(bliss_add_test module_name module_link)
    # set the name for the test target and executable
    set(test_target_name test-${module_name})
    # add all give files
    message(STATUS ${ARGN})
    add_executable(${test_target_name} ${ARGN})
    # add the googletest dependency
    add_dependencies(${test_target_name} googletest)
    # set binary output path for tests
    set_target_properties(${test_target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BINARY_OUTPUT_DIR})
    # link test with gtest_main
    target_link_libraries(${test_target_name} gtest_main)

    # all the extra libs
    target_link_libraries(${test_target_name} ${EXTRA_LIBS})

    # link to the tested module (but only if that module produces a linkable
    # library)
    if (module_link)
        target_link_libraries(${test_target_name} ${module_name})
    endif (module_link)
    
    # if code coverage is to be determined: link with gcov
    if (ENABLE_COVERAGE)
      target_link_libraries(${test_target_name} gcov)
    endif(ENABLE_COVERAGE)

    cotire(${test_target_name})

    # generate google test XML results, to be parsed by Jenkins
    # additionally to the CTest generated xml, the gtest xml is much more detailed (listing all single tests)
    add_test(NAME ${test_target_name} WORKING_DIRECTORY ${TEST_BINARY_OUTPUT_DIR} COMMAND ${test_target_name} ARGS "--gtest_output=xml:${TEST_XML_OUTPUT_DIR}/gtest/${test_target_name}.xml")
  endfunction(bliss_add_test)

  function(bliss_add_non_google_test module_names module_link)
      foreach(f ${ARGN})
        get_filename_component(test_name ${f} NAME_WE)
  
        # set the name for the test target and executable
        set(test_target_name ${test_name})
        # add all give files
        add_executable(${test_target_name} ${f})
        # set binary output path for tests
        set_target_properties(${test_target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BINARY_OUTPUT_DIR})
    
        # all the extra libs
        target_link_libraries(${test_target_name} ${EXTRA_LIBS})
    
        # link to the tested module (but only if that module produces a linkable
        # library)
        if (module_link)
            target_link_libraries(${test_target_name} ${module_names})
        endif (module_link)
        
        # if code coverage is to be determined: link with gcov
        if (ENABLE_COVERAGE)
          target_link_libraries(${test_target_name} gcov)
        endif(ENABLE_COVERAGE)
    
        cotire(${test_target_name})
    
        # generate google test XML results, to be parsed by Jenkins
        # additionally to the CTest generated xml, the gtest xml is much more detailed (listing all single tests)
        add_test(NAME ${test_target_name} WORKING_DIRECTORY ${TEST_BINARY_OUTPUT_DIR} COMMAND ${test_target_name})
      endforeach(f)
  
  endfunction(bliss_add_non_google_test)

  function(bliss_add_metaprog_test module_names module_link metadef)
      foreach(f ${ARGN})
        get_filename_component(test_name ${f} NAME_WE)
  
        # set the name for the test target and executable
        set(test_target_name ${test_name}_${metadef})
        # add all give files
        add_executable(${test_target_name} ${f})
        # set binary output path for tests
        set_target_properties(${test_target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BINARY_OUTPUT_DIR})
        SET_TARGET_PROPERTIES(${test_target_name} PROPERTIES COMPILE_FLAGS -D${metadef})
        
        # all the extra libs
        target_link_libraries(${test_target_name} ${EXTRA_LIBS})
    
        # link to the tested module (but only if that module produces a linkable
        # library)
        if (module_link)
            target_link_libraries(${test_target_name} ${module_names})
        endif (module_link)
        
        # if code coverage is to be determined: link with gcov
        if (ENABLE_COVERAGE)
          target_link_libraries(${test_target_name} gcov)
        endif(ENABLE_COVERAGE)
    
        cotire(${test_target_name})
    
        # generate google test XML results, to be parsed by Jenkins
        # additionally to the CTest generated xml, the gtest xml is much more detailed (listing all single tests)
        add_test(NAME ${test_target_name} WORKING_DIRECTORY ${TEST_BINARY_OUTPUT_DIR} COMMAND ${test_target_name})
      endforeach(f)
  
  endfunction(bliss_add_metaprog_test)
endif (ENABLE_TESTING)


###### Executable and Libraries
# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Output directory for applications" )

#link_directories( ${LIBRARY_OUTPUT_PATH} )



# Check if the user want to build test applications
CMAKE_DEPENDENT_OPTION(BUILD_TEST_APPLICATIONS "Inform whether test applications should be built" ON
                        "ENABLE_TESTING" OFF)
# Check if the user want to build sample applications
OPTION(BUILD_EXAMPLE_APPLICATIONS "Inform whether sample applications should be built" OFF)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config/config.hpp.in"
    "${PROJECT_BINARY_DIR}/config.hpp"
    @ONLY
    )
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config/logger_config.hpp.in"
    "${PROJECT_BINARY_DIR}/logger_config.hpp"
    @ONLY
    )


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/ext")
include_directories("${PROJECT_SOURCE_DIR}/ext/farmhash/src")
#include_directories("${PROJECT_SOURCE_DIR}/ext/sparsehash/src")




###### LIBRARY
#set(libsrcs CmdlineParser.cpp Communicator_I.cpp CommHandler_I.cpp Process.cpp PullCommHandler.cpp PushCommHandler.cpp RandomScheduler.cpp RoundRobinScheduler.cpp NullSinkAction.cpp DataBuffer.cpp MPIDataBuffer.cpp MPISendDataBuffer.cpp MPIRecvDataBuffer.cpp)
#set(libname runtime)
# Add rule to generate execution engine, as a shared library (.so)
#add_library("${libname}" SHARED ${libsrcs})
# Link opencv stuff to the execution engine

add_subdirectory(src/utils)

set(EXTRA_LIBS ${EXTRA_LIBS} bliss-utils)

add_subdirectory(src/taskrunner)
add_subdirectory(src/partition)
add_subdirectory(src/common)
add_subdirectory(src/concurrent)
add_subdirectory(src/io)
add_subdirectory(src/index)
add_subdirectory(src/iterators)
add_subdirectory(test)
add_subdirectory(utils)


# Install Settings


# Installer Settings


# Test Settings
