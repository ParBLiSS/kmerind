set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})
cotire(test_omp_patterns)

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

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP)
target_link_libraries(benchmark_concurrent_IO_MMAP ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_MMAP )


# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_ADVISE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_ADVISE)
target_link_libraries(benchmark_concurrent_IO_MMAP_ADVISE ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_MMAP_ADVISE )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_POPULATE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE)
target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_MMAP_POPULATE )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_ADVISE)
target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE )

# HUGETLB DOES NOT WORK WITH FILES.
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_HUGEPAGE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_HUGEPAGE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_HUGEPAGE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_HUGEPAGE ${EXTRA_LIBS})
#cotire(benchmark_concurrent_IO_MMAP_HUGEPAGE )
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_HUGEPAGE_ADVISE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE ${EXTRA_LIBS})
#cotire(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE )
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_HUGEPAGE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE ${EXTRA_LIBS})
#cotire(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE )
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_HUGEPAGE_ADVISE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE ${EXTRA_LIBS})
#cotire(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE )
#


# EXECUTABLES
add_executable(benchmark_concurrent_IO_FILELOADER benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FILELOADER PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER)
target_link_libraries(benchmark_concurrent_IO_FILELOADER ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FILELOADER )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FILELOADER_ATOMIC benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FILELOADER_ATOMIC PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER_ATOMIC)
target_link_libraries(benchmark_concurrent_IO_FILELOADER_ATOMIC ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FILELOADER_ATOMIC )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQ benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQ PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQ)
target_link_libraries(benchmark_concurrent_IO_FASTQ ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FASTQ )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIter benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIter PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter)
target_link_libraries(benchmark_concurrent_IO_FASTQIter ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FASTQIter )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIter2 benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIter2 PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter2)
target_link_libraries(benchmark_concurrent_IO_FASTQIter2 ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FASTQIter2 )

# EXECUTABLES
add_executable(benchmark_concurrent_IO_FASTQIterNoQual benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_FASTQIterNoQual PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIterNoQual)
target_link_libraries(benchmark_concurrent_IO_FASTQIterNoQual ${EXTRA_LIBS})
cotire(benchmark_concurrent_IO_FASTQIterNoQual )



add_executable(testKmerIndex_KMOLECULE TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMOLECULE PROPERTIES COMPILE_FLAGS -DKMOLECULEINDEX)
target_link_libraries(testKmerIndex_KMOLECULE ${EXTRA_LIBS})
cotire( testKmerIndex_KMOLECULE)

add_executable(testKmerIndex_KMER TestKmerIndexAndQuery.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_KMER PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_KMER ${EXTRA_LIBS})
cotire(testKmerIndex_KMER )

#add_executable(testKmerIndex_ST TestKmerIndexAndQuerySingleThread.cpp)
#SET_TARGET_PROPERTIES(testKmerIndex_ST PROPERTIES COMPILE_FLAGS -DKMERINDEX)
#target_link_libraries(testKmerIndex_ST ${EXTRA_LIBS})
#cotire(testKmerIndex_ST )
#
add_executable(testKmerIndex_ST2 TestKmerIndexAndQuerySingleThread_DistMap.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_ST2 PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_ST2 ${EXTRA_LIBS})
#cotire(testKmerIndex_ST2 )

add_executable(testKmerIndex_read TestKmerIndex_ParseOnly.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_read PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_read ${EXTRA_LIBS})
#cotire(testKmerIndex_ST2 )

add_executable(testKmerIndex_read_fasta TestKmerIndex_FASTAParseOnly.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_read_fasta PROPERTIES COMPILE_FLAGS -DKMERINDEX)
target_link_libraries(testKmerIndex_read_fasta ${EXTRA_LIBS})
#cotire(testKmerIndex_ST2 )


#add_executable(testKmerIndex_ST2 TestKmerIndexAndQuerySingleThread_DHashMap.cpp)
#SET_TARGET_PROPERTIES(testKmerIndex_ST2 PROPERTIES COMPILE_FLAGS -DKMERINDEX)
#target_link_libraries(testKmerIndex_ST2 ${EXTRA_LIBS})
#cotire(testKmerIndex_ST2 )
#
#add_executable(testKmerIndex_ST3 TestKmerIndexAndQuerySingleThread_DSortedMap.cpp)
#SET_TARGET_PROPERTIES(testKmerIndex_ST3 PROPERTIES COMPILE_FLAGS -DKMERINDEX)
#target_link_libraries(testKmerIndex_ST3 ${EXTRA_LIBS})
#cotire(testKmerIndex_ST3 )
#
#add_executable(testKmerIndex_ST4 TestKmerIndexAndQuerySingleThread_DMap.cpp)
#SET_TARGET_PROPERTIES(testKmerIndex_ST4 PROPERTIES COMPILE_FLAGS -DKMERINDEX)
#target_link_libraries(testKmerIndex_ST4 ${EXTRA_LIBS})
#cotire(testKmerIndex_ST4 )
#

add_executable(testKmerCountCorrectness TestKmerCountCorrectness.cpp)
target_link_libraries(testKmerCountCorrectness ${EXTRA_LIBS})

add_executable(testKmerCountCorrectnessInFasta TestKmerCountCorrectnessInFasta.cpp)
target_link_libraries(testKmerCountCorrectnessInFasta ${EXTRA_LIBS})

add_executable(test_de_brujin_graph_construction test_de_brujin_graph_construction.cpp)
target_link_libraries(test_de_brujin_graph_construction ${EXTRA_LIBS})

## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp) 
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})
cotire(profile_mem_usage )


add_executable(test_get_file_size test_get_file_size.cpp)
target_link_libraries(test_get_file_size ${EXTRA_LIBS})
cotire(test_get_file_size )


add_executable(hash_vs_sort hash_vs_sort.cpp)
target_link_libraries(hash_vs_sort ${EXTRA_LIBS})


add_executable(kmer_reverse kmer_reverse.cpp)
target_link_libraries(kmer_reverse ${EXTRA_LIBS})

if (Boost_FOUND)
	add_executable(cust_alloc test_custom_allocator.cpp)
	target_link_libraries(cust_alloc ${EXTRA_LIBS})
endif()

endif()

