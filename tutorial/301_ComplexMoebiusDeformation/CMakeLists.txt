cmake_minimum_required(VERSION 3.16)
project(301_ComplexMoebiusDeformation)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::glfw igl::opengl tutorials)
