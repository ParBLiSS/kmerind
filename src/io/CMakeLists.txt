project(bliss-io)
include("${CMAKE_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
    # get all mpi test files from ./test
    FILE(GLOB MPI_TEST_FILES test/mpi_test_*.cpp)
    bliss_add_mpi_test(${PROJECT_NAME} FALSE ${MPI_TEST_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    
    bliss_add_relacy_test(${PROJECT_NAME} FALSE test/relacy_messagebuffers.cpp BLISS_MUTEX)
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_messagebuffers.cpp BLISS_MUTEX)

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_messagebuffers2.cpp BLISS_MUTEX)

endif()



endif()
