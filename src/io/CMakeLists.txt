project(bliss-io)
include("${CMAKE_SOURCE_DIR}/Sanitizer.cmake")

if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB TEST_FILES test/test_*.cpp)
    bliss_add_test(${PROJECT_NAME} FALSE ${TEST_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
endif()

# load the testing:
if (IS_DIRECTORY ./test)
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE LOCKING test/check2_buffer.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE SPINLOCKING test/check2_buffer.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE LOCKFREE test/check2_buffer.cpp )
#    bliss_add_metaprog_test(${PROJECT_NAME} FALSE LOCKFREE_PTR test/check2_buffer.cpp )
endif()



endif()