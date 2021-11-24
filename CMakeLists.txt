

cmake_minimum_required(VERSION 3.1)
project(flat)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

add_subdirectory(lib/libigl)

file(GLOB SRCFILES src/*.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES} app/flat.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::png igl::opengl_glfw_imgui igl::triangle)
include_directories(include)

add_executable(procustes app/procustes.cpp)
target_link_libraries(procustes igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(net app/net.cpp src/net_param.cpp)
target_link_libraries(net igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(render_coords app/render_coords.cpp src/net_param.cpp)
target_link_libraries(render_coords igl::opengl_glfw igl::core)

