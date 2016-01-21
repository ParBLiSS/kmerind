set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(testCommLayer test_commlayer.cpp)
target_link_libraries(testCommLayer ${EXTRA_LIBS})
cotire(testCommLayer )

## EXECUTABLES
add_executable(testCommLayer_leaky test_commlayer_leakybarrier.cpp)
target_link_libraries(testCommLayer_leaky ${EXTRA_LIBS})
cotire(testCommLayer_leaky )

## EXECUTABLES
add_executable(stressTestCommLayer_sendonly stresstest_commlayer_sendonly.cpp)
target_link_libraries(stressTestCommLayer_sendonly ${EXTRA_LIBS})
cotire(stressTestCommLayer_sendonly )


## EXECUTABLES
add_executable(stressTestCommLayer_sendrecv stresstest_commlayer_sendrecv.cpp)
target_link_libraries(stressTestCommLayer_sendrecv ${EXTRA_LIBS})
cotire(stressTestCommLayer_sendrecv )


add_executable(testDistrMap test_distr_map.cpp)
target_link_libraries(testDistrMap ${EXTRA_LIBS})
cotire(testDistrMap )

endif()






if (BUILD_TEST_APPLICATIONS)


add_executable(testKmerIndex_KMOLECULE TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMOLECULE PROPERTIES COMPILE_FLAGS -DKMOLECULEINDEX)
target_link_libraries(testKmerIndex_KMOLECULE ${EXTRA_LIBS})
cotire( testKmerIndex_KMOLECULE)

add_executable(testKmerIndex_KMER TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMER PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_KMER ${EXTRA_LIBS})
cotire(testKmerIndex_KMER )

add_executable(testKmerCountCorrectness TestKmerCountCorrectness.cpp)
target_link_libraries(testKmerCountCorrectness ${EXTRA_LIBS})

#add_executable(testKmerCountCorrectnessInFasta TestKmerCountCorrectnessInFasta.cpp)
#target_link_libraries(testKmerCountCorrectnessInFasta ${EXTRA_LIBS})

add_executable(testKmerIndex_ST TestKmerIndexAndQuerySingleThread.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_ST PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_ST ${EXTRA_LIBS})
#cotire(testKmerIndex_ST )



endif()

