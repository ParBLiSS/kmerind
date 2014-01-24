cmake_minimum_required(VERSION 2.8)
project(pbil-common)

add_library(pbil-common
    dummy.cpp)

# target_link_libraries(pbil-common ...)

# load the testing:
if (IS_DIRECTORY ./test)
    message("found folder test")
else()
    message("noo test folder found")
endif()
