cmake_minimum_required(VERSION 3.14)
project(sac LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(sac ${SOURCES})

option(ENABLE_AVX2 "AVX2/FMA" OFF)
option(ENABLE_AVX512 "AVX512" OFF)

if(ENABLE_AVX512)
    target_compile_options(sac PRIVATE
        -mavx512f
        -mavx512dq
        -mavx512bw
        -mavx512vl
        -mfma
        -O3
    )
    target_compile_definitions(sac PRIVATE USE_AVX512=1)
elseif(ENABLE_AVX256)
    target_compile_options(sac PRIVATE
        -mavx2
        -mfma
        -O3
    )
    target_compile_definitions(sac PRIVATE USE_AVX256=1)
else()
    target_compile_options(sac PRIVATE
        -O3
    )
    target_compile_definitions(sac PRIVATE)
endif()

target_link_libraries(sac PRIVATE stdc++)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
