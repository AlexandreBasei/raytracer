# CMakeList.txt : projet CMake pour NewRaytracer, incluez la source et définissez
# la logique spécifique au projet ici.
#

# Ajoutez une source à l'exécutable de ce projet.
add_executable (NewRaytracer "RayTracer2.cpp" "external/rtw_stb_image.h" "external/stb_image.h" "perlin.h" "quad.h" "constant_medium.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET NewRaytracer PROPERTY CXX_STANDARD 20)
endif()

# Copy images directory to build output directory
add_custom_command(TARGET NewRaytracer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/images
    $<TARGET_FILE_DIR:NewRaytracer>/images
)

# TODO: Ajoutez des tests et installez des cibles si nécessaire.
