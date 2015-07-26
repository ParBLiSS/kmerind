
#add_subdirectory(testlib)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


## EXECUTABLES
add_executable(sys_probe sys_probe.cpp)
target_link_libraries(sys_probe ${EXTRA_LIBS})


add_executable(TextInspector text_inspector.cpp)
target_link_libraries(TextInspector ${EXTRA_LIBS})

