find_package(Python COMPONENTS Interpreter REQUIRED)
add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz"
  COMMAND ${Python_EXECUTABLE} -m build -o ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
  WORKING_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR})
  
add_custom_target(_python ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz")
