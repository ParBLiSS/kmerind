set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")



if (BUILD_TEST_APPLICATIONS)

#----------------
# COMMENTED OUT because the template parameters are not yet updated.
# also BenchmarkKmerIndex is serving the same purpose.
#add_executable(testKmerIndex_ST2 TestKmerIndexAndQuerySingleThread_DistMap.cpp)
#target_link_libraries(testKmerIndex_ST2 ${EXTRA_LIBS})

#add_executable(testKmerIndex_read TestKmerIndex_ParseOnly.cpp)
#target_link_libraries(testKmerIndex_read ${EXTRA_LIBS})

#add_executable(testKmerIndex_read_fasta TestKmerIndex_FASTAParseOnly.cpp)
#target_link_libraries(testKmerIndex_read_fasta ${EXTRA_LIBS})
#-----------------

add_executable(test_de_bruijn_graph_construction test_de_bruijn_graph_construction.cpp)
target_link_libraries(test_de_bruijn_graph_construction ${EXTRA_LIBS})

add_executable(test_get_file_size test_get_file_size.cpp)
target_link_libraries(test_get_file_size ${EXTRA_LIBS})

endif()

