cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(linefit)

find_package(PCL 1.8 REQUIRED)

include_directories(
  include
  ${PCL_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (${PROJECT_NAME}
  src/main.cpp
  src/linefit.cpp
)

target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES})
