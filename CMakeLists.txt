cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 20)

# List projects

# Libraries
project(simplex_method)
# Applications
project(lab_1)
project(lab_2)

# Enable Conan package manager
include(${CMAKE_CURRENT_SOURCE_DIR}/conan/conanbuildinfo.cmake)
set(CONAN_DISABLE_CHECK_COMPILER true)
conan_basic_setup()

# Library headers
set(LIB_SIMPLEX_METHOD_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/projects/simplex_method/include)

# Libraries

# [simplex_method]
file(GLOB LIB_SIMPLEX_METHOD_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/projects/simplex_method/src/*.cpp")
add_library(simplex_method STATIC ${LIB_SIMPLEX_METHOD_SOURCES})
target_include_directories(simplex_method PUBLIC ${LIB_SIMPLEX_METHOD_INCLUDE})

# Applications

function(create_lab_application lab_number)
    file(GLOB "LAB_${lab_number}_SOURCES" "${CMAKE_CURRENT_SOURCE_DIR}/projects/lab_${lab_number}/src/*.cpp")
    add_executable("lab_${lab_number}" "${LAB_${lab_number}_SOURCES}")
    target_include_directories("lab_${lab_number}" PUBLIC ${LIB_SIMPLEX_METHOD_SOURCES})
    target_link_libraries("lab_${lab_number}" simplex_method ${CONAN_LIBS})
endfunction()

# [lab_1]
CREATE_LAB_APPLICATION(1)
CREATE_LAB_APPLICATION(2)
