
#add_subdirectory(testlib)


## EXECUTABLES
add_executable(sys_probe sys_probe.cpp)
target_link_libraries(sys_probe ${EXTRA_LIBS})
#target_link_libraries(test testlib)


add_executable(TextInspector TextInspector.cpp)
target_link_libraries(TextInspector ${EXTRA_LIBS})

