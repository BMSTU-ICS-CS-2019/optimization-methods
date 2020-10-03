cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 20)

# List projects

# Libraries
project(simplex_method)
# Applications
project(lab_1)

# Enable Conan package manager
include(${CMAKE_CURRENT_SOURCE_DIR}/conan/conanbuildinfo.cmake)
conan_basic_setup()

# Library headers
set(LIB_SIMPLEX_METHOD_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/projects/simplex_method/include)

# Libraries

# [simplex_method]
file(GLOB LIB_SIMPLEX_METHOD_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/projects/simplex_method/src/*.cpp")
add_library(simplex_method STATIC ${LIB_SIMPLEX_METHOD_SOURCES})
target_include_directories(simplex_method PUBLIC ${LIB_SIMPLEX_METHOD_INCLUDE})

# Applications

# [lab_1]
file(GLOB LAB_1_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/projects/lab_1/src/*.cpp")
add_executable(lab_1 ${LAB_1_SOURCES})
target_include_directories(lab_1 PUBLIC ${LIB_SIMPLEX_METHOD_SOURCES})
target_link_libraries(lab_1 simplex_method ${CONAN_LIBS})
