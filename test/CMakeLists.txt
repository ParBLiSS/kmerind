
#add_subdirectory(testlib)


## EXECUTABLES
add_executable(test test.cpp)
target_link_libraries(test ${EXTRA_LIBS})


## EXECUTABLES
#add_executable(testIndex2 testIndex.cpp)
#target_link_libraries(testIndex2 ${EXTRA_LIBS})
#target_link_libraries(testIndex2 bliss-io)
#target_link_libraries(testIndex2 bliss-common)


## EXECUTABLES
add_executable(testTypes test_types.cpp)
target_link_libraries(testTypes ${EXTRA_LIBS})


## EXECUTABLES
add_executable(testThreads test_threads.cpp)
target_link_libraries(testThreads ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-io)
target_link_libraries(testThreads bliss-common)
#target_link_libraries(testThreads bliss-index)

## EXECUTABLES
add_executable(testThreadsNoQual test_threads-no_qual.cpp)
target_link_libraries(testThreadsNoQual ${EXTRA_LIBS})
#target_link_libraries(testThreadsNoQual bliss-io)
target_link_libraries(testThreadsNoQual bliss-common)
#target_link_libraries(testThreadsNoQual bliss-index)


## EXECUTABLES
add_executable(testIndexEl quicktest.cpp)
target_link_libraries(testIndexEl ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-io)
target_link_libraries(testIndexEl bliss-common)
#target_link_libraries(testThreads bliss-index)
