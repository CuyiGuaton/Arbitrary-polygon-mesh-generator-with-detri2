cmake_minimum_required(VERSION 3.0.0)
project(DelaunayPolyGenerator VERSION 0.1.0)

add_subdirectory(detri2)
add_executable(DelaunayPolyGenerator main.cpp polygon.cpp io.cpp triang.cpp)

target_link_libraries(DelaunayPolyGenerator detri2)



