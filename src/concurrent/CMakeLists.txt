project(bliss-concurrent)
include("${CMAKE_SOURCE_DIR}/cmake/Sanitizer.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/ExtraWarnings.cmake")


if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)

    bliss_add_relacy_test(${PROJECT_NAME} FALSE test/relacy_practice.cpp USE_RELACY  )

    bliss_add_relacy_test(${PROJECT_NAME} FALSE test/relacy_buffer.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE BLISS_NONE)
    bliss_add_relacy_test(${PROJECT_NAME} FALSE test/relacy_threadsafe_queue.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE )

    bliss_add_relacy_test(${PROJECT_NAME} FALSE test/relacy_objectpool.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE BLISS_NONE)

    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
    
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_buffer.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_buffer_performance.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE )
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_objectpool.cpp BLISS_MUTEX BLISS_SPINLOCK )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_objectpool_performance.cpp BLISS_MUTEX BLISS_SPINLOCK )

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_objectpool_lf.cpp BLISS_LOCKFREE )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_objectpool_performance_lf.cpp BLISS_LOCKFREE )

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_threadsafe_queue.cpp BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE test/check2_threadsafe_queue_performance.cpp BLISS_NONE BLISS_MUTEX BLISS_SPINLOCK BLISS_LOCKFREE )
    
endif()

endif()