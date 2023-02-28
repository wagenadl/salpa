set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles" ${CMAKE_MODULE_PATH})
find_package(Octave QUIET)

if (OCTAVE_MKOCTFILE)
    
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/salpa/private/salpamex.mex"
      COMMAND ${OCTAVE_MKOCTFILE} --mex -o "${CMAKE_CURRENT_BINARY_DIR}/salpa/private/salpamex.mex" salpamex.cpp LocalFit.cpp
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/salpa/private")
    
    configure_file("salpa/salpa.m" "salpa/salpa.m" COPYONLY)
    configure_file("salpa/private/getopt.m" "salpa/private/getopt.m" COPYONLY)
    configure_file("salpa/private/strtoks.m" "salpa/private/strtoks.m" COPYONLY)
    add_custom_target(_mex ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/salpa/private/salpamex.mex")

else()
  message(WARNING "Octave not found - not building octave module")
endif()
    
