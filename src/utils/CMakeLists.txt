project(bliss-utils)


## EXECUTABLES
add_library(${PROJECT_NAME} logging.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS})



# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} TRUE ${TEST_FILES})
endif()

