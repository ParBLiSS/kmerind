if (ENABLE_SANITIZER)
  # set flags for coverage test
  message(STATUS "Sanitizer enabled")
  # probably should be specific based on sanitizer style.  no-omit-frame-pointer is for address.
  add_definitions(-fsanitize=${SANITIZER_STYLE} -fPIE -fno-omit-frame-pointer)
  # and static-libtsan is for thread.
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=${SANITIZER_STYLE} -fno-omit-frame-pointer -pie -static-libtsan")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fsanitize=${SANITIZER_STYLE} -fno-omit-frame-pointer -pie -static-libtsan")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=${SANITIZER_STYLE} -fno-omit-frame-pointer -pie -static-libtsan")
  set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fsanitize=${SANITIZER_STYLE} -fno-omit-frame-pointer -pie -static-libtsan")
endif()
