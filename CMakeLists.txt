
# Pasos para compilar
# Generar CMakeLists.txt
# Crear carpeta build y acceder a ella
# Ejecutar 'cmake ..'
# Ejecutar 'make'

#cmake_minimum_required(VERSION 3.1...3.23)
#project( Quinto )

#find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

#find_package(CGAL REQUIRED COMPONENTS Qt5)
#add_definitions(-DCGAL_USE_BASIC_VIEWER)

#create_single_source_cgal_program( "cgal.cpp" )
#create_single_source_cgal_program( "voronoi.cpp" )

cmake_minimum_required(VERSION 3.14)
project(ActInt2)

# Find CGAL and its dependencies
find_package(CGAL REQUIRED)

# Create the executable
add_executable(ActInt2 probe2.cpp)

# Link standard CGAL libraries
target_link_libraries(ActInt2 PRIVATE CGAL::CGAL)