
#add_subdirectory(testlib)


## EXECUTABLES
add_executable(testIndex indexing.cpp)
target_link_libraries(testIndex ${EXTRA_LIBS})
#target_link_libraries(test testlib)
