
#add_subdirectory(testlib)


## EXECUTABLES
add_executable(sys_probe sys_probe.cpp)
target_link_libraries(sys_probe ${EXTRA_LIBS})
#target_link_libraries(test testlib)
