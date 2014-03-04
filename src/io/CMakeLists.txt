cmake_minimum_required(VERSION 2.8)
project(bliss-io)


## EXECUTABLES
add_library(${PROJECT_NAME} fastq_loader.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS})
#target_link_libraries(test testlib)



# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/*.cpp)
    bliss_add_test(${PROJECT_NAME} TRUE ${TEST_FILES})
endif()

