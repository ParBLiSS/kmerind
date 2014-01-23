
#add_subdirectory(testlib)


## EXECUTABLES
add_executable(test test.cpp)
target_link_libraries(test ${EXTRA_LIBS})
#target_link_libraries(test testlib)


## EXECUTABLES
add_executable(testFastQ testFastQ.cpp)
target_link_libraries(testFastQ ${EXTRA_LIBS})
#target_link_libraries(test testlib)
