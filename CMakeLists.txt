cmake_minimum_required(VERSION 3.10)
project(RayTracer VERSION 1.0 LANGUAGES CXX)

# Utiliser C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sources (ajoute/supprime les fichiers selon ton arborescence)
set(SOURCES
    ./RayTracer/RayTracer2.cpp
)

add_executable(RayTracer ${SOURCES})

# Include dirs si tu as des headers dans un sous-dossier
target_include_directories(RayTracer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Linker des libs externes si besoin
# target_link_libraries(RayTracer PRIVATE <lib1> <lib2>)