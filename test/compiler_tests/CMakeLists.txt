set(TEST_NAME bliss_compiler_tests)
include("${PROJECT_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/ExtraWarnings.cmake")

if(ENABLE_TESTING)


## EXECUTABLES
add_executable(gcc521-O3-loop gcc521-O3-loop.cpp)
target_link_libraries(gcc521-O3-loop ${EXTRA_LIBS})

add_executable(icc-autodeduce icc-autodeduce.cpp)
target_link_libraries(icc-autodeduce ${EXTRA_LIBS})

endif()



