project(bliss-index)
include("${CMAKE_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if (ENABLE_TESTING)
#add_subdirectory(testlib)
# TODO: add .cpp files (if we have some)
#add_library(${PROJECT_NAME}
#    kmer_index_functors.cpp)
#target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS})

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp test/benchmark_*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
    # get all mpi test files from ./test
    FILE(GLOB MPI_TEST_FILES test/mpi_test_*.cpp)
    bliss_add_mpi_test(${PROJECT_NAME} FALSE ${MPI_TEST_FILES})
endif()
endif()
