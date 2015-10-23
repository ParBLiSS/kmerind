set(TEST_NAME bliss-common)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${TEST_NAME} FALSE ${TEST_FILES})

    FILE(GLOB BENCHMARK_FILES test/benchmark_*.cpp)
    bliss_add_benchmark(${TEST_NAME} FALSE ${BENCHMARK_FILES})

endif()


endif()
