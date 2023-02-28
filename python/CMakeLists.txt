find_package(Python COMPONENTS Interpreter REQUIRED)
#find_package(PythonExtensions REQUIRED)
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()

add_library(salpapy MODULE src/salpapy.cpp src/LocalFit.cpp)
#python_extension_module(_salpapy)
install(TARGETS salpapy LIBRARY DESTINATION salpa)

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz"
    COMMAND ${Python_EXECUTABLE} -m build -o ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS salpapy)
  
add_custom_target(_python ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/salpa-1.0.tar.gz")
