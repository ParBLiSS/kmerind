set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})


endif()






if (BUILD_TEST_APPLICATIONS)


add_executable(testFASTQ_load BenchmarkFileLoader.cpp)
SET_TARGET_PROPERTIES(testFASTQ_load PROPERTIES COMPILE_FLAGS -DUSE_FASTQ_PARSER)
target_link_libraries(testFASTQ_load ${EXTRA_LIBS})

add_executable(testFASTA_load BenchmarkFileLoader.cpp)
SET_TARGET_PROPERTIES(testFASTA_load PROPERTIES COMPILE_FLAGS -DUSE_FASTA_PARSER)
target_link_libraries(testFASTA_load ${EXTRA_LIBS})


add_executable(testBuildIndex_fastq BenchmarkKmerIndexBuild.cpp)
SET_TARGET_PROPERTIES(testBuildIndex_fastq PROPERTIES COMPILE_FLAGS -DUSE_FASTQ_PARSER)
target_link_libraries(testBuildIndex_fastq ${EXTRA_LIBS})

add_executable(testBuildIndex_fasta BenchmarkKmerIndexBuild.cpp)
SET_TARGET_PROPERTIES(testBuildIndex_fasta PROPERTIES COMPILE_FLAGS -DUSE_FASTA_PARSER)
target_link_libraries(testBuildIndex_fasta ${EXTRA_LIBS})


if (USE_OPENMP)

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP)
target_link_libraries(benchmark_concurrent_IO_MMAP ${EXTRA_LIBS})


# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_ADVISE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_ADVISE)
target_link_libraries(benchmark_concurrent_IO_MMAP_ADVISE ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_POPULATE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE)
target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE ${EXTRA_LIBS})

# EXECUTABLES
add_executable(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE benchmark_concurrent_IO.cpp)
SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_ADVISE)
target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_ADVISE ${EXTRA_LIBS})

# HUGETLB DOES NOT WORK WITH FILES.
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_HUGEPAGE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_HUGEPAGE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_HUGEPAGE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_HUGEPAGE ${EXTRA_LIBS})
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_HUGEPAGE_ADVISE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_HUGEPAGE_ADVISE ${EXTRA_LIBS})
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_HUGEPAGE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE ${EXTRA_LIBS})
#
## EXECUTABLES
#add_executable(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE benchmark_concurrent_IO.cpp)
#SET_TARGET_PROPERTIES(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP_POPULATE_HUGEPAGE_ADVISE)
#target_link_libraries(benchmark_concurrent_IO_MMAP_POPULATE_HUGEPAGE_ADVISE ${EXTRA_LIBS})
#

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

endif(USE_OPENMP)


## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp) 
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})


add_executable(hash_vs_sort hash_vs_sort.cpp)
target_link_libraries(hash_vs_sort ${EXTRA_LIBS})

add_executable(chrono_vs_time chrono_vs_time.cpp)
target_link_libraries(chrono_vs_time ${EXTRA_LIBS} -lrt)


if (Boost_FOUND)
	add_executable(cust_alloc test_custom_allocator.cpp)
	target_link_libraries(cust_alloc ${EXTRA_LIBS})
endif()

endif()

