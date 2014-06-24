

## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})

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
add_executable(test_iterator_and_function_traits test_iterator_and_function_traits.cpp)
target_link_libraries(test_iterator_and_function_traits ${EXTRA_LIBS})


## EXECUTABLES
add_executable(testThreads test_threads.cpp)
target_link_libraries(testThreads ${EXTRA_LIBS})
target_link_libraries(testThreads bliss-common)



## EXECUTABLES
add_executable(profile_mem_usage profile_mem_usage.cpp)
target_link_libraries(profile_mem_usage ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-io)
target_link_libraries(profile_mem_usage bliss-common)
#target_link_libraries(testThreads bliss-index)
