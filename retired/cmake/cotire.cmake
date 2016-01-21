### tool to speed up compilation
if(CMAKE_VERSION VERSION_GREATER 2.8.3)
    OPTION(ENABLE_COTIRE "Enable using cotire to reduce compilation time (possibly)" OFF)
    if (ENABLE_COTIRE)
        set (CMAKE_MODULE_PATH ${EXT_PROJECTS_DIR}/cotire/CMake ${CMAKE_MODULE_PATH})
        include(cotire)
    ELSE()
        function(cotire module_names)
        endfunction(cotire)  
    ENDIF()
else()
    function(cotire module_names)
    endfunction(cotire)  
endif()

# call with
# cotire(target)