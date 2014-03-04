cmake_minimum_required(VERSION 2.8)
project(bliss-iterators)


# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
endif()

