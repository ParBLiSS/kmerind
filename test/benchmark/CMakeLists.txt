set(TEST_NAME bliss_test_apps)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(test_omp_patterns test_omp_patterns.cpp)
target_link_libraries(test_omp_patterns ${EXTRA_LIBS})


endif()

if (BL_BENCHMARK)

# available params:


# the following is tested as a group
# pDNA  (4, 5, 16)  -- affects any that uses LEX or XOR transform.
# pKmerStore  (SINGLE, CANONICAL, BIMOLECULE)
# pMAP  COUNT(SORTED, UNORDERED)  POS(SORTED, COMPACTVEC)

# test as group  SINGLE
# pDNA  (4, 5, 16)  -- affects any that uses LEX or XOR transform. 
# pDistTrans  (IDEN, XOR, LEX)   - only for pKmerStore=SINGLE and pMAP!=SORTED - test separately from above.
# pMAP count(unordered)   pos(compactvec)

# the following is tested as a group
# pK    (any...) - for scalability
# pKmerStore (SINGLE, CANONICAL, BIMOLECULE)

# pPARSER  (FASTA, FASTQ)        test FASTA separately
# pINDEX  (COUNT, POS, POSQUAL)  test POSQUAL separately.
# pMAP count(ORDERED)  POS(ORDERED UNORDERED VEC)-  test different backends separately.

# pDistHash (STD, IDEN, FARM, MURMUR) - NOT for pMAP=SORTED.  test separately
# pStoreHash (STD, IDEN, FARM, MURMUR) - NOT for pMAP=SORTED or pMAP=ORDERED.  test separately
# pCollective, pIrecv  ( turn on a2a or send-irecv based find)  test separately


# file naming:
#  testKmerIndex_${pPARSER}_a${pDNA}_k${pK}_${pKmerStore}_${pMAP}_${pINDEX}_XXXXXX
#  where XXXXXX are parameters are relevant for different Map types
#   these prefixed by 2 letters, and in this order of appearance: 
#       pDistTrans, pDistHash, pStoreHash, pCollective/pIrecv
#                           prefix          SORT      ORDERED     HASHED
#              pDistTrans:  dt              n         y           y           affects only singlestrand kmerstore
#               pDistHash:  DH              n         y           y
#              pStoreHash:  SH              n         n           y
#      pIrecv/pCollective:  CC              y         y           y

function(add_sortedmap_target file prefix parser dna k store map index)
  #  disttrans disthash storehash  ignored if passed in
  
      add_executable(${prefix}-${parser}-a${dna}-k${k}-${store}-SORTED-${index}-dtXXXX-dhYYYY-shZZZZ ${file})
      
      SET_TARGET_PROPERTIES(${prefix}-${parser}-a${dna}-k${k}-${store}-SORTED-${index}-dtXXXX-dhYYYY-shZZZZ
         PROPERTIES COMPILE_FLAGS 
         "-DpPARSER=${parser} -DpDNA=${dna} -DpK=${k} -DpKmerStore=${store} -DpMAP=SORTED -DpINDEX=${index}")
      
      target_link_libraries(${prefix}-${parser}-a${dna}-k${k}-${store}-SORTED-${index}-dtXXXX-dhYYYY-shZZZZ
       ${EXTRA_LIBS})


endfunction(add_sortedmap_target)

function(add_orderedmap_target file prefix parser dna k store map index disttrans disthash)
    # storehash is ignored if passed in

      add_executable(${prefix}-${parser}-a${dna}-k${k}-${store}-ORDERED-${index}-dt${disttrans}-dh${disthash}-shZZZZ ${file})
      
      SET_TARGET_PROPERTIES(${prefix}-${parser}-a${dna}-k${k}-${store}-ORDERED-${index}-dt${disttrans}-dh${disthash}-shZZZZ
         PROPERTIES COMPILE_FLAGS 
         "-DpPARSER=${parser} -DpDNA=${dna} -DpK=${k} -DpKmerStore=${store} -DpMAP=ORDERED -DpINDEX=${index} -DpDistTrans=${disttrans} -DpDistHash=${disthash}")
      
      target_link_libraries(${prefix}-${parser}-a${dna}-k${k}-${store}-ORDERED-${index}-dt${disttrans}-dh${disthash}-shZZZZ
       ${EXTRA_LIBS})

endfunction(add_orderedmap_target)

function(add_hashmap_target file prefix parser dna k store map index disttrans disthash storehash)

      add_executable(${prefix}-${parser}-a${dna}-k${k}-${store}-${map}-${index}-dt${disttrans}-dh${disthash}-sh${storehash} ${file})
      
      SET_TARGET_PROPERTIES(${prefix}-${parser}-a${dna}-k${k}-${store}-${map}-${index}-dt${disttrans}-dh${disthash}-sh${storehash}
         PROPERTIES COMPILE_FLAGS 
         "-DpPARSER=${parser} -DpDNA=${dna} -DpK=${k} -DpKmerStore=${store} -DpMAP=${map} -DpINDEX=${index} -DpDistTrans=${disttrans} -DpDistHash=${disthash} -DpStoreHash=${storehash}")
      
      target_link_libraries(${prefix}-${parser}-a${dna}-k${k}-${store}-${map}-${index}-dt${disttrans}-dh${disthash}-sh${storehash}
       ${EXTRA_LIBS})


endfunction(add_hashmap_target)

#foreach(store SINGLE CANONICAL BIMOLECULE)
#	add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} HASHEDVEC POS IDEN FARM FARM)
#endforeach(store)

#================= 36 + 12 targets
# vary alphabet, which affects kmer transformations and has impact on type of backend containers.
# test all 3 kmerstores at the default hash and transforms. with different types of maps.

foreach(dna 4 5 16)
  
  # base set. - vary alphabet and kmer store, then check count and pos for sorted vs hashed.
  foreach(store SINGLE CANONICAL BIMOLECULE)
  
    # count maps.  note SORTED PATH ignores hash but uses transformation
    add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} SORTED COUNT IDEN FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} UNORDERED COUNT IDEN FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} DENSEHASH COUNT IDEN FARM FARM)
    
    # position maps.  note SORTED PATH ignores hash but uses transformation
    add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} SORTED POS IDEN FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} UNORDERED POS IDEN FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 ${store} DENSEHASH POS IDEN FARM FARM)
  endforeach(store)
  
  # then add for SINGLE strand store, different distribution transforms for COUNT and POS
  foreach(disttrans XOR LEX)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 SINGLE UNORDERED COUNT ${disttrans} FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 SINGLE DENSEHASH COUNT ${disttrans} FARM FARM)
    
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 SINGLE UNORDERED POS ${disttrans} FARM FARM)
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ ${dna} 31 SINGLE DENSEHASH POS ${disttrans} FARM FARM)
  endforeach(disttrans)

endforeach(dna)

#=================  108 targets, varying K
#vary K, fix Alphabet, and map types.  see effect on count and position (map vs multimap), and different kmerstores
# do it for fastA and fastQ

foreach(store SINGLE CANONICAL BIMOLECULE)

	  foreach(format FASTA FASTQ)
	    foreach(k 15 21 63 95) 
	        # count maps.  note SORTED PATH ignores hash but uses transformation
	        add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex ${format} 4 ${k} ${store} SORTED COUNT IDEN FARM FARM)
	     	add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex ${format} 4 ${k} ${store} DENSEHASH COUNT IDEN FARM FARM)
	        
	        # position maps.  note SORTED PATH ignores hash but uses transformation
	        add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex ${format} 4 ${k} ${store} SORTED POS IDEN FARM FARM)
	    	add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex ${format} 4 ${k} ${store} DENSEHASH POS IDEN FARM FARM)
	    endforeach(k)
	  endforeach(format)

	        # count maps.  note SORTED PATH ignores hash but uses transformation
	        add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTA 4 31 ${store} SORTED COUNT IDEN FARM FARM)
	     	add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTA 4 31 ${store} DENSEHASH COUNT IDEN FARM FARM)
	        
	        # position maps.  note SORTED PATH ignores hash but uses transformation
	        add_sortedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTA 4 31 ${store} SORTED POS IDEN FARM FARM)
	    	add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTA 4 31 ${store} DENSEHASH POS IDEN FARM FARM)

endforeach(store)


#==================  1 target  quality map.  Single to assess quality effect.
# pos quality maps.. 
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE DENSEHASH POSQUAL IDEN FARM FARM)

    
#================== 8 targets - slow backends, or potentially no advantage
# store model changes the collision characteristics, so study these...
foreach(store SINGLE BIMOLECULE) 
  add_orderedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} ORDERED COUNT IDEN FARM FARM)
  
  add_orderedmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} ORDERED POS IDEN FARM FARM)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} VEC POS IDEN FARM FARM)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} COMPACTVEC POS IDEN FARM FARM)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 ${store} HASHEDVEC POS IDEN FARM FARM)
endforeach(store)

#=====================  6  targets
# vary distribution hash method.  use SINGLE to reduce collision due to lex_less.

foreach(hash IDEN STD MURMUR)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE DENSEHASH COUNT IDEN ${hash} FARM)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE DENSEHASH POS IDEN ${hash} FARM)
endforeach(hash)  

#=====================  18  targets
# vary storage hash method. use SINGLE to reduce collision due to lex_less.

foreach(hash IDEN STD MURMUR)
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE UNORDERED COUNT IDEN FARM ${hash})
  add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE DENSEHASH COUNT IDEN FARM ${hash})
  foreach(map VEC COMPACTVEC UNORDERED DENSEHASH HASHEDVEC) 
    add_hashmap_target(BenchmarkKmerIndex.cpp testKmerIndex FASTQ 4 31 SINGLE ${map} POS IDEN FARM ${hash})
  endforeach(map)
endforeach(hash)  


#================== 8 targets,  pos only (need the data volume)  - potentially slow hence their own separate tests.
# collective find tests.  want high collision (imbalance) to see how all2all handles.
# no count.  no VEC - very similar to compactvec. fix a, k.  don't worry about xor.  compare BIMOLE vs SINGLE.

foreach(store BIMOLECULE SINGLE)

  #add_executable(testKmerIndex-COMPQ-a4-k31-${store}-VEC-POS-dtIDEN-dhFARM-shFARM BenchmarkKmerIndex.cpp)
  #SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-VEC-POS-dtIDEN-dhFARM-shFARM
  #   PROPERTIES COMPILE_FLAGS 
  #   "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=VEC -DpINDEX=POS -DpDistTrans=IDEN -DpDistHash=FARM -DpStoreHash=FARM -DpCompare")
  #target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-VEC-POS-dtIDEN-dhFARM-shFARM
  # ${EXTRA_LIBS})
  
#  add_executable(testKmerIndex-COMPQ-a4-k31-${store}-UNORDERED-POS-dtIDEN-dhFARM-shFARM BenchmarkKmerIndex.cpp)
#  SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-UNORDERED-POS-dtIDEN-dhFARM-shFARM
#     PROPERTIES COMPILE_FLAGS 
#     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=UNORDERED -DpINDEX=POS -DpDistTrans=IDEN -DpDistHash=FARM -DpStoreHash=FARM -DpCompare")
#  target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-UNORDERED-POS-dtIDEN-dhFARM-shFARM
#   ${EXTRA_LIBS})
  
  #add_executable(testKmerIndex-COMPQ-a4-k31-${store}-ORDERED-POS-dtIDEN-dhFARM-shZZZZ BenchmarkKmerIndex.cpp)
  #SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-ORDERED-POS-dtIDEN-dhFARM-shZZZZ
  #   PROPERTIES COMPILE_FLAGS 
  #   "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=ORDERED -DpINDEX=POS -DpDistTrans=IDEN -DpDistHash=FARM -DpCompare")
  #target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-ORDERED-POS-dtIDEN-dhFARM-shZZZZ
  # ${EXTRA_LIBS})
  
  add_executable(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-POS-dtIDEN-dhFARM-shFARM BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-POS-dtIDEN-dhFARM-shFARM
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=DENSEHASH -DpINDEX=POS -DpDistTrans=IDEN -DpDistHash=FARM -DpStoreHash=FARM -DpCompare")
  target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-POS-dtIDEN-dhFARM-shFARM
   ${EXTRA_LIBS})

  add_executable(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-COUNT-dtXXXX-dhYYYY-shZZZZ BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-COUNT-dtXXXX-dhYYYY-shZZZZ
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=DENSEHASH -DpINDEX=COUNT -DpCompare")
  target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-DENSEHASH-COUNT-dtXXXX-dhYYYY-shZZZZ
   ${EXTRA_LIBS})


  add_executable(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-POS-dtIDEN-dhFARM-shFARM BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-POS-dtIDEN-dhFARM-shFARM
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=SORTED -DpINDEX=POS -DpDistTrans=IDEN -DpDistHash=FARM -DpStoreHash=FARM -DpCompare")
  target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-POS-dtIDEN-dhFARM-shFARM
   ${EXTRA_LIBS})

  add_executable(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-COUNT-dtXXXX-dhYYYY-shZZZZ BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-COUNT-dtXXXX-dhYYYY-shZZZZ
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=${store} -DpMAP=SORTED -DpINDEX=COUNT -DpCompare")
  target_link_libraries(testKmerIndex-COMPQ-a4-k31-${store}-SORTED-COUNT-dtXXXX-dhYYYY-shZZZZ
   ${EXTRA_LIBS})
 
endforeach(store)

#==================


endif()

if (BUILD_TEST_APPLICATIONS)

  add_executable(testKmerIndex-fastq BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-fastq
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTQ -DpDNA=4 -DpK=31 -DpKmerStore=BIMOLECULE -DpMAP=DENSEHASH -DpINDEX=POSQUAL -DpDistTrans=IDEN -DpDistHash=MURMUR -DpStoreHash=FARM -DpCompare")
  target_link_libraries(testKmerIndex-fastq
   ${EXTRA_LIBS})

  add_executable(testKmerIndex-fasta BenchmarkKmerIndex.cpp)
  SET_TARGET_PROPERTIES(testKmerIndex-fasta
     PROPERTIES COMPILE_FLAGS 
     "-DpPARSER=FASTA -DpDNA=5 -DpK=21 -DpKmerStore=CANONICAL -DpMAP=SORTED -DpINDEX=COUNT -DpCompare")
  target_link_libraries(testKmerIndex-fasta
   ${EXTRA_LIBS})


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


# EXECUTABLES
include_directories("${EXT_PROJECTS_DIR}/tommyds")
add_executable(benchmark_hashtables BenchmarkHashTables.cpp)
target_link_libraries(benchmark_hashtables ${EXTRA_LIBS})
