
#add_subdirectory(testlib)


## EXECUTABLES
#add_executable(test1 test.cpp)
#target_link_libraries(test1 ${EXTRA_LIBS})

add_executable(TextInspector TextInspector.cpp)
target_link_libraries(TextInspector ${EXTRA_LIBS})


# EXECUTABLES
add_executable(testCIO_MMAP testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_MMAP PROPERTIES COMPILE_FLAGS -DTEST_OP_MMAP)
target_link_libraries(testCIO_MMAP ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FILELOADER testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FILELOADER PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER)
target_link_libraries(testCIO_FILELOADER ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FILELOADER_ATOMIC testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FILELOADER_ATOMIC PROPERTIES COMPILE_FLAGS -DTEST_OP_FILELOADER_ATOMIC)
target_link_libraries(testCIO_FILELOADER_ATOMIC ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FASTQ testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FASTQ PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQ)
target_link_libraries(testCIO_FASTQ ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FASTQIter testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FASTQIter PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter)
target_link_libraries(testCIO_FASTQIter ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FASTQIter2 testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FASTQIter2 PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIter2)
target_link_libraries(testCIO_FASTQIter2 ${EXTRA_LIBS})

# EXECUTABLES
add_executable(testCIO_FASTQIterNoQual testCIO.cpp)
SET_TARGET_PROPERTIES(testCIO_FASTQIterNoQual PROPERTIES COMPILE_FLAGS -DTEST_OP_FASTQIterNoQual)
target_link_libraries(testCIO_FASTQIterNoQual ${EXTRA_LIBS})



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
target_link_libraries(testThreads bliss-common)

### EXECUTABLES
#add_executable(testThreadsNoQual test_threads-no_qual.cpp)
#target_link_libraries(testThreadsNoQual ${EXTRA_LIBS})
##target_link_libraries(testThreadsNoQual bliss-io)
#target_link_libraries(testThreadsNoQual bliss-common)
##target_link_libraries(testThreadsNoQual bliss-index)
#
#
### EXECUTABLES
#add_executable(testThreadsNoMPI test_threads_no_mpi.cpp)
#target_link_libraries(testThreadsNoMPI ${EXTRA_LIBS})
##target_link_libraries(testThreadsNoMPI bliss-io)
#target_link_libraries(testThreadsNoMPI bliss-common)
##target_link_libraries(testThreadsNoMPI bliss-index)
#
#
### EXECUTABLES
#add_executable(testThreadsNoQualNoMPI test_threads_no_mpi_no_qual.cpp)
#target_link_libraries(testThreadsNoQualNoMPI ${EXTRA_LIBS})
##target_link_libraries(testThreadsNoQualNoMPI bliss-io)
#target_link_libraries(testThreadsNoQualNoMPI bliss-common)
##target_link_libraries(testThreadsNoQualNoMPI bliss-index)

## EXECUTABLES
add_executable(testIndexEl quicktest.cpp)
target_link_libraries(testIndexEl ${EXTRA_LIBS})
#target_link_libraries(testThreads bliss-io)
target_link_libraries(testIndexEl bliss-common)
#target_link_libraries(testThreads bliss-index)
