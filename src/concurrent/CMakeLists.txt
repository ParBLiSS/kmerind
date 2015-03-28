project(bliss-concurrent)
include("${bliss_SOURCE_DIR}/cmake/Sanitizer.cmake")


if (ENABLE_TESTING)
# load the testing:
if (IS_DIRECTORY ./test)
    # get all files from ./test
    FILE(GLOB CHECK_FILES test/check_*.cpp)
    
    bliss_add_non_google_test(${PROJECT_NAME} FALSE ${CHECK_FILES})
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_buffer.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_buffer.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_buffer.cpp )

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_buffer_performance.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_buffer_performance.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_buffer_performance.cpp )

    
    
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_objectpool.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_objectpool.cpp )

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_objectpool_performance.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_objectpool_performance.cpp )


    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_objectpool_lf.cpp )

    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_objectpool_performance_lf.cpp )


#    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_threadsafe_queue.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_threadsafe_queue.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_threadsafe_queue.cpp )


    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_NONE test/check2_threadsafe_queue_performance.cpp )
#    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_MUTEX test/check2_threadsafe_queue_performance.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_SPINLOCK test/check2_threadsafe_queue_performance.cpp )
    bliss_add_metaprog_test(${PROJECT_NAME} FALSE BLISS_LOCKFREE test/check2_threadsafe_queue_performance.cpp )
    
endif()

endif()
