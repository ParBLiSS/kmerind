set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})


endif()

if (BL_BENCHMARK)

#================= 72 targets
# vary alphabet, which affects kmer transformations as well.  this in turn has impact on type of backend containers.
foreach(dna 4 5 16)
  foreach(trans IDEN LEX XOR PRELEX)
  
    # count maps.  note SORTED PATH ignores hash but uses transformation
    foreach(map SORTED UNORDERED) 
      add_executable(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_COUNT BenchmarkKmerIndex.cpp)
      SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_COUNT
         PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=${dna} -DpK=31 -DpTRANS=${trans} -DpHASH=FARM -DpMAP=${map} -DpINDEX=COUNT")
      target_link_libraries(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_COUNT ${EXTRA_LIBS})
    endforeach(map)
    
    # position maps.  note SORTED PATH ignores hash but uses transformation
    foreach(map SORTED VEC COMPACTVEC UNORDERED) 
      add_executable(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_POS BenchmarkKmerIndex.cpp)
      SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_POS
         PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=${dna} -DpK=31 -DpTRANS=${trans} -DpHASH=FARM -DpMAP=${map} -DpINDEX=POS")
      target_link_libraries(testKmerIndex_FASTQ_A${dna}_K31_${trans}_FARM_${map}_POS ${EXTRA_LIBS})
    endforeach(map)
    
  endforeach(trans)

endforeach(dna)

#=================  16 targets
#vary K, fix Alphabet, hash, transform, and map types.  do for count and position (map vs multimap)
foreach(k 15 21 63 95)  # 31 already compiled.

    # count maps.  note SORTED PATH ignores hash but uses transformation
    foreach(map SORTED UNORDERED) 
      add_executable(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_COUNT BenchmarkKmerIndex.cpp)
      SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_COUNT
         PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=${k} -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=${map} -DpINDEX=COUNT")
      target_link_libraries(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_COUNT ${EXTRA_LIBS})
    endforeach(map)
    
    # position maps.  note SORTED PATH ignores hash but uses transformation
    foreach(map SORTED COMPACTVEC) 
      add_executable(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_POS BenchmarkKmerIndex.cpp)
      SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_POS
         PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=${k} -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=${map} -DpINDEX=POS")
      target_link_libraries(testKmerIndex_FASTQ_A4_K${k}_PRELEX_FARM_${map}_POS ${EXTRA_LIBS})
    endforeach(map)
  
endforeach(k)

#=====================  12 targets
# vary hash method.  FARM already compiled.  check ordered and unordered

foreach(hash IDEN STD MURMUR)

  # count maps
  foreach(map UNORDERED) 
    add_executable(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_COUNT BenchmarkKmerIndex.cpp)
    SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_COUNT
       PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=${hash} -DpMAP=${map} -DpINDEX=COUNT")
    target_link_libraries(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_COUNT ${EXTRA_LIBS})
  endforeach(map)
  
  # position maps
  foreach(map VEC COMPACTVEC UNORDERED) 
    add_executable(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_POS BenchmarkKmerIndex.cpp)
    SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_POS
       PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=${hash} -DpMAP=${map} -DpINDEX=POS")
    target_link_libraries(testKmerIndex_FASTQ_A4_K31_PRELEX_${hash}_${map}_POS ${EXTRA_LIBS})
  endforeach(map)
  
endforeach(hash)

#=================  2 target
# FASTA 
add_executable(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_UNORDERED_COUNT BenchmarkKmerIndex.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_UNORDERED_COUNT 
  PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTA -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=UNORDERED -DpINDEX=COUNT")
target_link_libraries(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_UNORDERED_COUNT ${EXTRA_LIBS})

add_executable(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_COMPACTVEC_POS BenchmarkKmerIndex.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_COMPACTVEC_POS 
  PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTA -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=COMPACTVEC -DpINDEX=POS")
target_link_libraries(testKmerIndex_FASTA_A4_K31_PRELEX_FARM_COMPACTVEC_POS ${EXTRA_LIBS})

#================== 2 targets
# ORDERED
add_executable(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_COUNT BenchmarkKmerIndex.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_COUNT
   PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=ORDERED -DpINDEX=COUNT")
target_link_libraries(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_COUNT ${EXTRA_LIBS})

add_executable(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_POS BenchmarkKmerIndex.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_POS
   PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=ORDERED -DpINDEX=POS")
target_link_libraries(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_ORDERED_POS ${EXTRA_LIBS})

#==================  1 target
# pos quality maps.. 
add_executable(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_COMPACTVEC_POSQUAL BenchmarkKmerIndex.cpp)
SET_TARGET_PROPERTIES(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_COMPACTVEC_POSQUAL 
  PROPERTIES COMPILE_FLAGS "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpTRANS=PRELEX -DpHASH=FARM -DpMAP=COMPACTVEC -DpINDEX=POSQUAL")
target_link_libraries(testKmerIndex_FASTQ_A4_K31_PRELEX_FARM_COMPACTVEC_POSQUAL ${EXTRA_LIBS})

#==================


endif()

if (BUILD_TEST_APPLICATIONS)

add_executable(testFASTQ_load BenchmarkFileLoader.cpp)
SET_TARGET_PROPERTIES(testFASTQ_load PROPERTIES COMPILE_FLAGS -DUSE_FASTQ_PARSER)
target_link_libraries(testFASTQ_load ${EXTRA_LIBS})

add_executable(testFASTA_load BenchmarkFileLoader.cpp)
SET_TARGET_PROPERTIES(testFASTA_load PROPERTIES COMPILE_FLAGS -DUSE_FASTA_PARSER)
target_link_libraries(testFASTA_load ${EXTRA_LIBS})



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

