###### GCCFilt support, using my fixed version of gccfilter.
OPTION(ENABLE_GCCFILTER "Use GCCFilter for prettifying compiler output.  Requires Perl" OFF)

CMAKE_DEPENDENT_OPTION(GCCFILTER_IDE_COMPAT "If compiling from an IDE, turn off colorization for compatibility" OFF "ENABLE_GCCFILTER" OFF)

if (ENABLE_GCCFILTER)

if (GCCFILTER_IDE_COMPAT)
    set(GCCFILTER_COLORIZE "" CACHE INTERNAL "gccfilter flag for colorization" FORCE)
else (GCCFILTER_IDE_COMPAT)
    set(GCCFILTER_COLORIZE "-c" CACHE INTERNAL "gccfilter flag for colorization" FORCE)
endif(GCCFILTER_IDE_COMPAT)

SET(GCCFILTER_VERBOSITY 0 CACHE STRING "G++ STL message condensation level:  -1, 0, 1, 2, 3, 4.  0 is most terse, 4 is most verbose, -1 is unmodified.")

if("${GCCFILTER_VERBOSITY}" STREQUAL "4")
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_SOURCE_DIR}/utils/gccfilter/gccfilter ${GCCFILTER_COLORIZE}")
elseif("${GCCFILTER_VERBOSITY}" STREQUAL "3")
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_SOURCE_DIR}/utils/gccfilter/gccfilter ${GCCFILTER_COLORIZE} -r -w")
elseif("${GCCFILTER_VERBOSITY}" STREQUAL "2")
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_SOURCE_DIR}/utils/gccfilter/gccfilter ${GCCFILTER_COLORIZE} -r -w -i")
elseif("${GCCFILTER_VERBOSITY}" STREQUAL "1")
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_SOURCE_DIR}/utils/gccfilter/gccfilter ${GCCFILTER_COLORIZE} -a -i")
elseif("${GCCFILTER_VERBOSITY}" STREQUAL "0")
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_SOURCE_DIR}/utils/gccfilter/gccfilter ${GCCFILTER_COLORIZE} -a -w -i")    
endif()

endif(ENABLE_GCCFILTER)
 