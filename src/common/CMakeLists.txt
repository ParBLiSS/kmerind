cmake_minimum_required(VERSION 2.8)
project(bliss-common)

# TODO: add .cpp files (if we have some)
#add_library(${PROJECT_NAME}
#    )

# target_link_libraries(${PROJECT_NAME} ...)


# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/*.cpp)
    bliss_add_test(${PROJECT_NAME} ${TEST_FILES})
endif()
