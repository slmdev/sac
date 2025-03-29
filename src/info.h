#pragma once

#define SAC_VERSION "0.7.17"

#define TOSTRING_HELPER(x) #x
#define TOSTRING(x) TOSTRING_HELPER(x)

#ifdef __clang__
#  define COMPILER                                                             \
    "clang " TOSTRING(__clang_major__) "." TOSTRING(__clang_minor__            \
    ) "." TOSTRING(__clang_patchlevel__)
#elif defined(__GNUC__)
#  define COMPILER                                                             \
    "gcc" TOSTRING(__GNUC__) "." TOSTRING(__GNUC_MINOR__                       \
    ) "." TOSTRING(__GNUC_PATCHLEVEL__)
#else
#  define COMPILER "Unknown"
#endif

#ifdef __x86_64__
#  define ARCHITECTURE "x86_64"
#elifdef __i386__
#  define ARCHITECTURE "x86_32"
#else
#  define ARCHITECTURE "Unknown"
#endif

#if defined(USE_AVX256)
#  define AVX2_STATE "ON"
#else
#  define AVX2_STATE "OFF"
#endif