cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 20)

# List projects
project(lab_1)

# Enable Conan package manager
include(${CMAKE_CURRENT_SOURCE_DIR}/conan/conanbuildinfo.cmake)
conan_basic_setup()

# [lab_1]
file(GLOB LAB_1_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/projects/lab_1/src/*.cpp")
add_executable(lab_1 ${LAB_1_SOURCES})
target_link_libraries(lab_1 ${CONAN_LIBS})
