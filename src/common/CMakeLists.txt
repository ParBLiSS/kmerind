cmake_minimum_required(VERSION 2.8)
project(pbil-common)

add_library(${PROJECT_NAME}
    dummy.cpp)

# target_link_libraries(${PROJECT_NAME} ...)


# load the testing:
if (IS_DIRECTORY ./test)
    # include the test subdirectory
    #add_subdirectory(test)
    FILE(GLOB TEST_FILES test/*.cpp)
    pbil_add_test(${PROJECT_NAME} ${TEST_FILES})
endif()
