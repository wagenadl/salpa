# CMakeLists.txt - Part of SALPA, (C) Daniel Wagenaar 2023

cmake_minimum_required(VERSION 3.22)

project(SALPA VERSION 1.0.0)
set(CMAKE_PROJECT_HOMEPAGE_URL "http://github.com/wagenadl/salpa")

######################################################################
# Specify the C++ standard [must be done _before_ adding executable!]
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(salpa WIN32 MACOSX_BUNDLE src/salpa.cpp src/LocalFit.cpp)

add_subdirectory("docs")
add_subdirectory("python")
add_subdirectory("matlab")
