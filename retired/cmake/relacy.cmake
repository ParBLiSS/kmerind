
#### turn on/off relacy for thread testing.
CMAKE_DEPENDENT_OPTION(ENABLE_RELACY "Enable Thread testing with Relacy" OFF
                        "ENABLE_TESTING;NOT ENABLE_SANITIZER" OFF)
IF(ENABLE_RELACY)
    include_directories("${EXT_PROJECTS_DIR}/relacy_2_4")
    include_directories("${EXT_PROJECTS_DIR}/cq_extracted")
ENDIF()
                        
IF(ENABLE_RELACY)  
    function(bliss_add_relacy_test module_names module_link filename)
      get_filename_component(test_name ${filename} NAME_WE)
      
      foreach(metadef ${ARGN})
  
        # set the name for the test target and executable
        set(test_target_name ${test_name}_${metadef})
        # add all give files
        add_executable(${test_target_name} ${filename})
        # set binary output path for tests
        set_target_properties(${test_target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BINARY_OUTPUT_DIR})
        
        #get_target_property(targetCompDefs ${test_target_name} COMPILE_DEFINITIONS )
        set_target_properties(${test_target_name} PROPERTIES COMPILE_DEFINITIONS "${metadef} -DUSE_RELACY")

        # all the extra libs
        target_link_libraries(${test_target_name} ${EXTRA_LIBS})
    
        # link to the tested module (but only if that module produces a linkable
        # library)
        if (module_link)
            target_link_libraries(${test_target_name} ${module_names})
        endif (module_link)
        
        # if code coverage is to be determined: link with gcov
        if (ENABLE_COVERAGE)
          target_link_libraries(${test_target_name} gcov)
        endif(ENABLE_COVERAGE)
    
        # generate google test XML results, to be parsed by Jenkins
        # additionally to the CTest generated xml, the gtest xml is much more detailed (listing all single tests)
        add_test(NAME ${test_target_name} WORKING_DIRECTORY ${TEST_BINARY_OUTPUT_DIR} COMMAND ${test_target_name})
      endforeach(metadef)
  
  endfunction(bliss_add_relacy_test)
ELSE()
    function(bliss_add_relacy_test module_names module_link filename)
  endfunction(bliss_add_relacy_test)

ENDIF()                       