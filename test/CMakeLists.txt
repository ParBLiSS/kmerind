project(bliss_test_apps)

include("${bliss_SOURCE_DIR}/cmake/Sanitizer.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})

## EXECUTABLES
add_executable(testCommLayer test_commlayer.cpp)
target_link_libraries(testCommLayer ${EXTRA_LIBS})
target_link_libraries(testCommLayer bliss-common)

### EXECUTABLES
#add_executable(testCommLayer_THREADSHARED test_commlayer.cpp)
#target_link_libraries(testCommLayer_THREADSHARED ${EXTRA_LIBS})
#target_link_libraries(testCommLayer_THREADSHARED bliss-common)


## EXECUTABLES
add_executable(stressTestCommLayer_sendonly stresstest_commlayer_sendonly.cpp)
target_link_libraries(stressTestCommLayer_sendonly ${EXTRA_LIBS})
target_link_libraries(stressTestCommLayer_sendonly bliss-common)

### EXECUTABLES
#add_executable(stressTestCommLayer_sendonly_THREADSHARED stresstest_commlayer_sendonly.cpp)
#target_link_libraries(stressTestCommLayer_sendonly_THREADSHARED ${EXTRA_LIBS})
#target_link_libraries(stressTestCommLayer_sendonly_THREADSHARED bliss-common)

## EXECUTABLES
add_executable(stressTestCommLayer_sendrecv stresstest_commlayer_sendrecv.cpp)
target_link_libraries(stressTestCommLayer_sendrecv ${EXTRA_LIBS})
target_link_libraries(stressTestCommLayer_sendrecv bliss-common)

### EXECUTABLES
#add_executable(stressTestCommLayer_sendrecv_THREADSHARED stresstest_commlayer_sendrecv.cpp)
#target_link_libraries(stressTestCommLayer_sendrecv_THREADSHARED ${EXTRA_LIBS})
#target_link_libraries(stressTestCommLayer_sendrecv_THREADSHARED bliss-common)


add_executable(testDistrMap test_distr_map.cpp)
target_link_libraries(testDistrMap ${EXTRA_LIBS})
target_link_libraries(testDistrMap bliss-common)

endif()






if (BUILD_TEST_APPLICATIONS)

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP)
target_link_libraries(benchmark_concurrent_IO_MMAP ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FILELOADER benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FILELOADER PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER)
target_link_libraries(benchmark_concurrent_IO_FILELOADER ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FILELOADER_ATOMIC benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FILELOADER_ATOMIC PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER_ATOMIC)
target_link_libraries(benchmark_concurrent_IO_FILELOADER_ATOMIC ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQ benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQ PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQ)
target_link_libraries(benchmark_concurrent_IO_FASTQ ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIter benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIter PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter)
target_link_libraries(benchmark_concurrent_IO_FASTQIter ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIter2 benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIter2 PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter2)
target_link_libraries(benchmark_concurrent_IO_FASTQIter2 ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIterNoQual benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIterNoQual PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIterNoQual)
target_link_libraries(benchmark_concurrent_IO_FASTQIterNoQual ${EXTRA_LIBS})



## EXECUTABLES
#add_executable(testIndex2 testIndex.cpp)
#target_link_libraries(testIndex2 ${EXTRA_LIBS})
#target_link_libraries(testIndex2 bliss-io)
#target_link_libraries(testIndex2 bliss-common)



## EXECUTABLES
#add_executable(testThreads test_threads.cpp)
#target_link_libraries(testThreads ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-common)

add_executable(testKmerIndex_KMOLECULE TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMOLECULE PROPERTIES COMPILE_FLAGS -DKMOLECULEINDEX)
target_link_libraries(testKmerIndex_KMOLECULE ${EXTRA_LIBS})
target_link_libraries(testKmerIndex_KMOLECULE bliss-common)

add_executable(testKmerIndex_KMER TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMER PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_KMER ${EXTRA_LIBS})
target_link_libraries(testKmerIndex_KMER bliss-common)



#add_executable(testKmerIndexNew TestKmerIndexAndQueryNew.cpp)
#target_link_libraries(testKmerIndexNew ${EXTRA_LIBS})
#target_link_libraries(testKmerIndexNew bliss-common)


add_executable(testKmerCountCorrectness TestKmerCountCorrectness.cpp)
target_link_libraries(testKmerCountCorrectness ${EXTRA_LIBS})
target_link_libraries(testKmerCountCorrectness bliss-common)

add_executable(testKmerCountCorrectnessInFasta TestKmerCountCorrectnessInFasta.cpp)
target_link_libraries(testKmerCountCorrectnessInFasta ${EXTRA_LIBS})
target_link_libraries(testKmerCountCorrectnessInFasta bliss-common)

## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp)
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-io)
target_link_libraries(profile_mem_usage bliss-common)
#target_link_libraries(testThreads bliss-index)



add_executable(quicktest quicktest.cpp)
target_link_libraries(quicktest ${EXTRA_LIBS})

endif()

