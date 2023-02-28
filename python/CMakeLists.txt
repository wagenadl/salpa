find_package(Python COMPONENTS Interpreter QUIET)

if (Python_EXECUTABLE)
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz"
      COMMAND ${Python_EXECUTABLE} -m build -o ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
      WORKING_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR})
      
    add_custom_target(_python ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz")
else()
  message(WARNING "Python not found - not building python module")
endif()    

