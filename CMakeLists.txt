

cmake_minimum_required(VERSION 3.1)
project(flat)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

find_package(OpenMP REQUIRED)

add_subdirectory(lib/libigl)

add_subdirectory(lib/json)
include_directories(lib/json/include)

file(GLOB SRCFILES src/*.cpp)
file(GLOB PARAFASHIONFILES src/parafashion/*.cpp)
file(GLOB TRACINGFILES src/field_tracing/*.cpp)

add_executable(${PROJECT_NAME} ${SRCFILES} app/flat.cpp src/param/param_utils.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::png igl::opengl_glfw_imgui igl::triangle)
include_directories(include)

add_executable(procustes app/procustes_app.cpp)
target_link_libraries(procustes igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(net app/net.cpp src/net_param.cpp)
target_link_libraries(net igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(render_coords app/render_coords.cpp src/net_param.cpp)
target_link_libraries(render_coords igl::opengl_glfw igl::core)

add_executable(triangle_energies app/triangle_energies.cpp src/net_param.cpp src/param/param_utils.cpp)
target_link_libraries(triangle_energies igl::opengl_glfw igl::core igl::opengl_glfw_imgui)

add_executable(mccartney_energies app/mccartney_energies.cpp src/net_param.cpp src/param/param_utils.cpp src/param/bary_optimizer.cpp src/param/auto_select.cpp)
target_link_libraries(mccartney_energies igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::png OpenMP::OpenMP_CXX igl::triangle)

add_executable(energy_check app/energy_check.cpp src/param/param_utils.cpp)
target_link_libraries(energy_check igl::core)

add_executable(compute_bnd app/compute_bnd.cpp src/param/param_utils.cpp)
target_link_libraries(compute_bnd igl::opengl_glfw igl::core igl::opengl_glfw_imgui)

add_executable(cut app/cut.cpp)
target_link_libraries(cut igl::core)

file(GLOB PARAMSRC src/param/*.cpp)

add_executable(measure_alignment app/measure_alignment.cpp ${PARAMSRC})
target_link_libraries(measure_alignment igl::core)

add_executable(multi_param app/multi_param.cpp ${PARAMSRC} src/param/metrics.cpp)
target_link_libraries(multi_param igl::core igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::triangle nlohmann_json)


add_executable(multiple_poses app/multiple_poses.cpp ${PARAMSRC})
target_link_libraries(multiple_poses igl::core igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::triangle nlohmann_json)

add_executable(param app/param.cpp ${PARAMSRC})
target_link_libraries(param igl::opengl_glfw igl::core igl::opengl_glfw_imgui)
target_include_directories(param PUBLIC include/param)

add_executable(make_uv app/make_uv.cpp)
target_link_libraries(make_uv igl::core)