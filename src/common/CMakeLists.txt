project(bliss-common)
include("${CMAKE_SOURCE_DIR}/cmake/Sanitizer.cmake")

# TODO: add .cpp files (if we have some)
add_library(${PROJECT_NAME}
    alphabets.cpp)

# target_link_libraries(${PROJECT_NAME} ...)

if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp test/benchmark_*.cpp)
    bliss_add_test(${PROJECT_NAME} TRUE ${TEST_FILES})
endif()
endif()