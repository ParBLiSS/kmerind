include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


## EXECUTABLES
add_library(bliss-utils logging.cpp)
target_link_libraries(bliss-utils ${EXTRA_LIBS})

if (ENABLE_TESTING)

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} TRUE ${TEST_FILES})
    # get all mpi test files from ./test
    FILE(GLOB MPI_TEST_FILES test/mpi_test_*.cpp)
    bliss_add_mpi_test(${PROJECT_NAME} TRUE ${MPI_TEST_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
endif()

endif()
