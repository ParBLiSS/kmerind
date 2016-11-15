set(TEST_NAME bliss-io)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${TEST_NAME} FALSE ${TEST_FILES})
    # get all mpi test files from ./test
    FILE(GLOB MPI_TEST_FILES test/mpi_test_*.cpp)
    bliss_add_mpi_test(${TEST_NAME} FALSE ${MPI_TEST_FILES})
if (BL_BENCHMARK)
    FILE(GLOB BENCHMARK_FILES test/benchmark_*.cpp)
    bliss_add_benchmark(${TEST_NAME} FALSE ${BENCHMARK_FILES})

    FILE(GLOB MPI_BENCHMARK_FILES test/mpi_benchmark_*.cpp)
    bliss_add_mpi_benchmark(${TEST_NAME} FALSE ${MPI_BENCHMARK_FILES})
endif()
    
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${TEST_NAME} FALSE ${CHECK_FILES})
endif()



endif()
