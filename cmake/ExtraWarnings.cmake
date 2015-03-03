

if(ENABLE_EXTRA_COMPILER_WARNINGS)
    add_definitions(${EXTRA_WARNING_FLAGS})
endif()


if(ENABLE_TYPE_CONVERSION_WARNINGS)
    add_definitions(${TYPE_CONVERSION_WARNING_FLAGS})
    
endif()

OPTION(ENABLE_SUGGESTION_WARNINGS "Enable compiler to generate suggestions for polymorphic type analysis?" OFF)
if(ENABLE_SUGGESTION_WARNINGS)
    add_definitions(${SUGGESTION_WARNING_FLAGS})
endif()

