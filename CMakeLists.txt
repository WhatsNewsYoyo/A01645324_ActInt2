cmake_minimum_required(VERSION 3.14)
project(ActInt2)

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_package(CGAL COMPONENTS Qt6 REQUIRED)

#Exe name
add_executable(ActInt2 main.cpp)

target_compile_definitions(ActInt2 PRIVATE CGAL_USE_BASIC_VIEWER)

#QT6 Libraries
target_link_libraries(ActInt2 PRIVATE
    CGAL::CGAL
    CGAL::CGAL_Qt6
)
