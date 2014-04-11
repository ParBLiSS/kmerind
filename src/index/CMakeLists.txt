project(bliss-index)


#add_subdirectory(testlib)
# TODO: add .cpp files (if we have some)
add_library(${PROJECT_NAME}
    kmer_index_functors.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS})


## EXECUTABLES
add_executable(testIndex indexing.cpp)
target_link_libraries(testIndex ${EXTRA_LIBS})
#target_link_libraries(test testlib)
