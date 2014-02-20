
#add_subdirectory(testlib)


## EXECUTABLES
add_library(bliss_io fastq_loader.cpp)
target_link_libraries(bliss_io ${EXTRA_LIBS})
#target_link_libraries(test testlib)

