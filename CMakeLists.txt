cmake_minimum_required(VERSION 3.10)
project(OBC_SC)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Packages
find_package(Eigen3 REQUIRED) 
include(FetchContent)
FetchContent_Declare(
        lime
        GIT_REPOSITORY https://github.com/usqcd-software/c-lime.git
        GIT_TAG        master  # Ou un tag spécifique comme v1.3.2
)
FetchContent_MakeAvailable(lime)



# Cible pointant vers le dossier src/
add_executable(test_gauge_geo src/bin/test_gauge_geo.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/su3/utils.cpp
    src/observables/observables.cpp
)
target_link_libraries(test_gauge_geo PRIVATE Eigen3::Eigen)

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(test_gauge_geo PRIVATE
            -O3
            -xHost       
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_gauge_geo PRIVATE
            -O3
            -ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()

add_executable(test_heatbath src/bin/test_heatbath.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/su3/utils.cpp
    src/observables/observables.cpp
    src/heatbath/heatbath.cpp
)
target_link_libraries(test_heatbath PRIVATE Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(test_heatbath PRIVATE
            -O3
            -xHost     
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_heatbath PRIVATE
            -O3
            -ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()

add_executable(test_ecmc src/bin/test_ecmc.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/su3/utils.cpp
    src/observables/observables.cpp
    src/ecmc/ecmc.cpp
    src/io/params.h
)
target_link_libraries(test_ecmc PRIVATE Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(test_ecmc PRIVATE
            -O3
            -xHost     
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_ecmc PRIVATE
            -O3
            -ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()

add_executable(gauge_ecmc src/bin/gauge_ecmc.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/su3/utils.cpp
    src/observables/observables.cpp
    src/ecmc/ecmc.cpp
    src/io/params.h
    src/io/io.cpp
    src/io/ildg.cpp
)
target_link_libraries(gauge_ecmc PRIVATE Eigen3::Eigen lime)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(gauge_ecmc PRIVATE
            -O3
            -xHost    
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_ecmc PRIVATE
            -O3
            -ffast-math
            -fopenmp
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()
