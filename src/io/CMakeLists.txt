project(bliss-io)
include("${CMAKE_SOURCE_DIR}/Sanitizer.cmake")

if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
endif()
endif()