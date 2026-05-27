cmake_minimum_required(VERSION 3.10)
project(OBC_SC)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED) 

# Cible pointant vers le dossier src/
add_executable(test_gauge_geo src/bin/test_gauge_geo.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/su3/utils.cpp
    src/observables/observables.cpp
)
target_link_libraries(test_gauge_geo PRIVATE Eigen3::Eigen)
