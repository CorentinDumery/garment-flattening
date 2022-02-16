

cmake_minimum_required(VERSION 3.1)
project(woven_param)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

find_package(OpenMP REQUIRED)

add_subdirectory(lib/libigl)
add_subdirectory(lib/implot)
add_subdirectory(lib/json)
include_directories(include)

file(GLOB PARAMSRC src/param/*.cpp)

add_library(${PROJECT_NAME} ${PARAMSRC})
target_include_directories(${PROJECT_NAME} PUBLIC include include/param)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core nlohmann_json)

add_executable(woven_viz app/woven_viz.cpp)
target_link_libraries(woven_viz woven_param igl::opengl_glfw igl::opengl_glfw_imgui igl::png)

add_executable(reflec_param app/reflec_param.cpp)
target_link_libraries(reflec_param woven_param igl::opengl_glfw igl::opengl_glfw_imgui igl::triangle)

add_executable(multiple_poses app/multiple_poses.cpp)
target_link_libraries(multiple_poses woven_param igl::opengl_glfw igl::opengl_glfw_imgui igl::png)

add_executable(procrustes_test app/procrustes_test.cpp)
target_link_libraries(procrustes_test woven_param igl::opengl_glfw_imgui)

add_executable(triangle_energies app/triangle_energies.cpp)
target_link_libraries(triangle_energies woven_param igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(measure_alignment app/measure_alignment.cpp)
target_link_libraries(measure_alignment woven_param)
