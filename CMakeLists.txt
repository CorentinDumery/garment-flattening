

cmake_minimum_required(VERSION 3.1)
project(woven_param)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

find_package(OpenMP REQUIRED)

add_subdirectory(lib/libigl)
add_subdirectory(lib/implot)
add_subdirectory(lib/json)
include_directories(include)

file(GLOB SRCFILES src/*.cpp)
file(GLOB PARAMSRC src/param/*.cpp)

add_library(${PROJECT_NAME} ${PARAMSRC})
target_include_directories(${PROJECT_NAME} PUBLIC include include/param)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core nlohmann_json)

add_executable(flat app/flat.cpp)
target_link_libraries(flat woven_param igl::core igl::opengl_glfw igl::png igl::opengl_glfw_imgui igl::triangle implot)

add_executable(procrustes_test appnew/procrustes_test.cpp)
target_link_libraries(procrustes_test woven_param igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(triangle_energies appnew/triangle_energies.cpp)
target_link_libraries(triangle_energies woven_param igl::opengl_glfw igl::core igl::opengl_glfw_imgui)

add_executable(mccartney_energies app/mccartney_energies.cpp src/param/param_utils.cpp src/param/bary_optimizer.cpp src/param/auto_select.cpp)
target_link_libraries(mccartney_energies igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::png OpenMP::OpenMP_CXX igl::triangle)

add_executable(measure_alignment appnew/measure_alignment.cpp)
target_link_libraries(measure_alignment woven_param)

add_executable(multi_param app/multi_param.cpp ${PARAMSRC} src/param/metrics.cpp)
target_link_libraries(multi_param igl::core igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::triangle nlohmann_json)

add_executable(multiple_poses appnew/multiple_poses.cpp ${PARAMSRC})
target_link_libraries(multiple_poses igl::core igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::png nlohmann_json)

add_executable(woven_viz appnew/woven_viz.cpp)
target_link_libraries(woven_viz woven_param igl::core igl::opengl_glfw igl::png igl::core igl::opengl_glfw_imgui nlohmann_json)