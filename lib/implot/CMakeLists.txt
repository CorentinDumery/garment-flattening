cmake_minimum_required(VERSION 3.1)

project(implot LANGUAGES CXX)

file(GLOB SRCIMPLOT src/*.cpp)

add_library(implot ${SRCIMPLOT})
target_link_libraries(implot PUBLIC igl::opengl_glfw_imgui)
target_include_directories(implot PUBLIC include)