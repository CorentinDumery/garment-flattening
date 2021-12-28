

cmake_minimum_required(VERSION 3.1)
project(flat)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

find_package(OpenMP REQUIRED)

add_subdirectory(lib/libigl)
#add_subdirectory(lib/vcglib)

file(GLOB SRCFILES src/*.cpp)
file(GLOB PARAFASHIONFILES src/parafashion/*.cpp)
file(GLOB TRACINGFILES src/field_tracing/*.cpp)

file(GLOB VCGSRC lib/vcglib/wrap/ply/plylib.cpp)
file(GLOB OLDIGLSRC lib/oldigl/igl/n_polyvector.cpp # oh my
                    lib/oldigl/igl/principal_curvature.cpp
                    lib/oldigl/igl/edge_topology.cpp
                    lib/oldigl/igl/per_face_normals.cpp
                    lib/oldigl/igl/per_vertex_normals.cpp
                    lib/oldigl/igl/vertex_triangle_adjacency.cpp
                    lib/oldigl/igl/slice.cpp
                    lib/oldigl/igl/internal_angles.cpp
                    lib/oldigl/igl/polyroots.cpp
                    lib/oldigl/igl/local_basis.cpp
                    lib/oldigl/igl/adjacency_list.cpp)

add_executable(${PROJECT_NAME} ${SRCFILES} app/flat.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::png igl::opengl_glfw_imgui igl::triangle)
include_directories(include)

add_executable(procustes app/procustes_app.cpp)
target_link_libraries(procustes igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(net app/net.cpp src/net_param.cpp)
target_link_libraries(net igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

add_executable(render_coords app/render_coords.cpp src/net_param.cpp)
target_link_libraries(render_coords igl::opengl_glfw igl::core)

add_executable(triangle_energies app/triangle_energies.cpp src/net_param.cpp src/procustes.cpp)
target_link_libraries(triangle_energies igl::opengl_glfw igl::core igl::opengl_glfw_imgui)

add_executable(mccartney_energies app/mccartney_energies.cpp src/net_param.cpp src/procustes.cpp src/bary_optimizer.cpp)
target_link_libraries(mccartney_energies igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::png OpenMP::OpenMP_CXX)

add_executable(energy_check app/energy_check.cpp src/procustes.cpp)
target_link_libraries(energy_check igl::core)

#add_executable(markables app/markables.cpp ${VCGSRC} ${OLDIGLSRC} ${PARAFASHIONFILES} ${TRACINGFILES})
#target_link_libraries(markables igl::opengl_glfw igl::core igl::opengl_glfw_imgui igl::png)
#target_include_directories(markables 
#PUBLIC lib/vcglib/    
#        lib/oldigl/    
#        lib/field_tracing/tracing/)

#add_executable(crosses app/crosses.cpp)
#target_link_libraries(crosses igl::opengl_glfw igl::core igl::opengl_glfw_imgui)
