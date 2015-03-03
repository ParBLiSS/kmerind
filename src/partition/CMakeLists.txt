project(bliss-partition)
include("${CMAKE_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
endif()

endif()